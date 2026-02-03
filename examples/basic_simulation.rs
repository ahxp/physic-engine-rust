//! Basic physics simulation example
//!
//! This example demonstrates a ball falling onto a floor under gravity.

use rustphy::prelude::*;

fn main() {
    println!("RustPhy - Basic Simulation Example");
    println!("===================================\n");

    // Create physics world with default settings
    let mut world = World::default();
    world.set_gravity(Vec3::new(0.0, -9.81, 0.0));

    // Create a static floor
    let floor = world.create_body(RigidBodyDesc::fixed().with_position(Vec3::ZERO));
    world.attach_collider(floor, Shape::cuboid(Vec3::new(10.0, 0.5, 10.0)), 1.0);
    println!("Created floor at Y=0 (top surface at Y=0.5)");

    // Create a dynamic ball
    let ball = world.create_body(
        RigidBodyDesc::dynamic()
            .with_position(Vec3::new(0.0, 5.0, 0.0))
            .with_mass(1.0),
    );
    world.attach_collider(ball, Shape::sphere(0.5), 1.0);
    println!("Created ball at Y=5.0 (radius=0.5)\n");

    // Simulation parameters
    let dt = 1.0 / 60.0;
    let total_time = 3.0;
    let steps = (total_time / dt) as usize;

    println!("Simulating {} seconds ({} steps at {}Hz)...\n", total_time, steps, 1.0 / dt);

    // Run simulation
    for i in 0..steps {
        world.step(dt);

        // Print position every 30 frames (0.5 seconds)
        if i % 30 == 0 {
            let pos = world.body_position(ball);
            let vel = world.body(ball).map(|b| b.linear_velocity).unwrap_or(Vec3::ZERO);
            println!(
                "t={:.2}s: position=({:.3}, {:.3}, {:.3}), velocity=({:.3}, {:.3}, {:.3})",
                i as f32 * dt,
                pos.x, pos.y, pos.z,
                vel.x, vel.y, vel.z
            );
        }
    }

    let final_pos = world.body_position(ball);
    println!("\nFinal ball position: ({:.3}, {:.3}, {:.3})", final_pos.x, final_pos.y, final_pos.z);
    println!("Expected resting position: ~(0, 1.0, 0) (floor top at 0.5 + ball radius 0.5)");
}
