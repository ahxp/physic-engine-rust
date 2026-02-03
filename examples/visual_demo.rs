//! Visual demo of the physics engine using macroquad
//!
//! Run with: cargo run --example visual_demo

use macroquad::prelude::*;
use rustphy::collision::BodyHandle;
use rustphy::dynamics::RigidBodyDesc;
use rustphy::geometry::Shape;
use rustphy::math::Vec3 as PhysVec3;
use rustphy::math::Quat as PhysQuat;
use rustphy::World;

// Window settings
const WINDOW_WIDTH: f32 = 1000.0;
const WINDOW_HEIGHT: f32 = 700.0;

// Physics to screen conversion (pixels per meter)
const SCALE: f32 = 40.0;

// Convert physics Y (up) to screen Y (down)
fn physics_to_screen(pos: PhysVec3) -> (f32, f32) {
    let x = WINDOW_WIDTH / 2.0 + pos.x * SCALE;
    let y = WINDOW_HEIGHT - 80.0 - pos.y * SCALE;
    (x, y)
}

fn window_conf() -> Conf {
    Conf {
        window_title: "RustPhy - Visual Demo".to_owned(),
        window_width: WINDOW_WIDTH as i32,
        window_height: WINDOW_HEIGHT as i32,
        ..Default::default()
    }
}

// Static platform info for drawing
struct Platform {
    handle: BodyHandle,
    half_extents: PhysVec3,
    color: Color,
}

#[macroquad::main(window_conf)]
async fn main() {
    // Create physics world
    let mut world = World::default();
    world.set_gravity(PhysVec3::new(0.0, -15.0, 0.0));

    // Store platforms for drawing
    let mut platforms: Vec<Platform> = Vec::new();

    // === Create ground floor ===
    let floor = world.create_body(
        RigidBodyDesc::fixed()
            .with_position(PhysVec3::new(0.0, 0.0, 0.0)),
    );
    let floor_extents = PhysVec3::new(12.0, 0.3, 1.0);
    world.attach_collider(floor, Shape::cuboid(floor_extents), 1.0);
    platforms.push(Platform { handle: floor, half_extents: floor_extents, color: DARKGRAY });

    // === Left ramp (tilted platform) ===
    let left_ramp = world.create_body(
        RigidBodyDesc::fixed()
            .with_position(PhysVec3::new(-6.0, 4.0, 0.0))
            .with_rotation(PhysQuat::from_axis_angle(PhysVec3::Z, 0.3)), // ~17 degrees
    );
    let ramp_extents = PhysVec3::new(3.0, 0.2, 1.0);
    world.attach_collider(left_ramp, Shape::cuboid(ramp_extents), 1.0);
    platforms.push(Platform { handle: left_ramp, half_extents: ramp_extents, color: Color::from_rgba(100, 100, 150, 255) });

    // === Right ramp (tilted opposite) ===
    let right_ramp = world.create_body(
        RigidBodyDesc::fixed()
            .with_position(PhysVec3::new(6.0, 4.0, 0.0))
            .with_rotation(PhysQuat::from_axis_angle(PhysVec3::Z, -0.3)),
    );
    world.attach_collider(right_ramp, Shape::cuboid(ramp_extents), 1.0);
    platforms.push(Platform { handle: right_ramp, half_extents: ramp_extents, color: Color::from_rgba(100, 100, 150, 255) });

    // === Middle platform ===
    let mid_platform = world.create_body(
        RigidBodyDesc::fixed()
            .with_position(PhysVec3::new(0.0, 7.0, 0.0)),
    );
    let mid_extents = PhysVec3::new(2.5, 0.2, 1.0);
    world.attach_collider(mid_platform, Shape::cuboid(mid_extents), 1.0);
    platforms.push(Platform { handle: mid_platform, half_extents: mid_extents, color: Color::from_rgba(80, 120, 80, 255) });

    // === Upper left platform ===
    let upper_left = world.create_body(
        RigidBodyDesc::fixed()
            .with_position(PhysVec3::new(-5.0, 10.0, 0.0)),
    );
    let small_extents = PhysVec3::new(1.5, 0.15, 1.0);
    world.attach_collider(upper_left, Shape::cuboid(small_extents), 1.0);
    platforms.push(Platform { handle: upper_left, half_extents: small_extents, color: Color::from_rgba(150, 80, 80, 255) });

    // === Upper right platform ===
    let upper_right = world.create_body(
        RigidBodyDesc::fixed()
            .with_position(PhysVec3::new(5.0, 10.0, 0.0)),
    );
    world.attach_collider(upper_right, Shape::cuboid(small_extents), 1.0);
    platforms.push(Platform { handle: upper_right, half_extents: small_extents, color: Color::from_rgba(150, 80, 80, 255) });

    // === Funnel walls ===
    let funnel_left = world.create_body(
        RigidBodyDesc::fixed()
            .with_position(PhysVec3::new(-3.0, 2.0, 0.0))
            .with_rotation(PhysQuat::from_axis_angle(PhysVec3::Z, -0.5)),
    );
    let funnel_extents = PhysVec3::new(1.5, 0.15, 1.0);
    world.attach_collider(funnel_left, Shape::cuboid(funnel_extents), 1.0);
    platforms.push(Platform { handle: funnel_left, half_extents: funnel_extents, color: Color::from_rgba(120, 100, 80, 255) });

    let funnel_right = world.create_body(
        RigidBodyDesc::fixed()
            .with_position(PhysVec3::new(3.0, 2.0, 0.0))
            .with_rotation(PhysQuat::from_axis_angle(PhysVec3::Z, 0.5)),
    );
    world.attach_collider(funnel_right, Shape::cuboid(funnel_extents), 1.0);
    platforms.push(Platform { handle: funnel_right, half_extents: funnel_extents, color: Color::from_rgba(120, 100, 80, 255) });

    // === Create dynamic balls ===
    let mut balls: Vec<(BodyHandle, f32, Color)> = Vec::new();

    // Large balls
    for i in 0..3 {
        let ball = world.create_body(
            RigidBodyDesc::dynamic()
                .with_position(PhysVec3::new(
                    (i as f32 - 1.0) * 2.0,
                    12.0 + i as f32 * 1.0,
                    0.0,
                ))
                .with_mass(2.0),
        );
        world.attach_collider(ball, Shape::sphere(0.6), 1.0);
        balls.push((ball, 0.6, RED));
    }

    // Medium balls
    for i in 0..4 {
        let ball = world.create_body(
            RigidBodyDesc::dynamic()
                .with_position(PhysVec3::new(
                    (i as f32 - 1.5) * 1.5,
                    15.0 + i as f32 * 0.8,
                    0.0,
                ))
                .with_mass(1.0),
        );
        world.attach_collider(ball, Shape::sphere(0.4), 1.0);
        balls.push((ball, 0.4, ORANGE));
    }

    // Small balls
    for i in 0..6 {
        let ball = world.create_body(
            RigidBodyDesc::dynamic()
                .with_position(PhysVec3::new(
                    (i as f32 - 2.5) * 1.0,
                    18.0 + i as f32 * 0.5,
                    0.0,
                ))
                .with_mass(0.5),
        );
        world.attach_collider(ball, Shape::sphere(0.25), 1.0);
        balls.push((ball, 0.25, YELLOW));
    }

    // === Create dynamic boxes ===
    let mut boxes: Vec<(BodyHandle, PhysVec3, Color)> = Vec::new();

    // Various sized boxes
    let box_configs = [
        (PhysVec3::new(-4.0, 14.0, 0.0), PhysVec3::new(0.5, 0.5, 0.5), BLUE, 2.0),
        (PhysVec3::new(4.0, 14.0, 0.0), PhysVec3::new(0.5, 0.5, 0.5), BLUE, 2.0),
        (PhysVec3::new(0.0, 16.0, 0.0), PhysVec3::new(0.7, 0.3, 0.5), PURPLE, 1.5),
        (PhysVec3::new(-2.0, 17.0, 0.0), PhysVec3::new(0.3, 0.6, 0.5), SKYBLUE, 1.0),
        (PhysVec3::new(2.0, 17.0, 0.0), PhysVec3::new(0.4, 0.4, 0.5), VIOLET, 1.2),
    ];

    for (pos, extents, color, mass) in box_configs {
        let box_body = world.create_body(
            RigidBodyDesc::dynamic()
                .with_position(pos)
                .with_mass(mass),
        );
        world.attach_collider(box_body, Shape::cuboid(extents), 1.0);
        boxes.push((box_body, extents, color));
    }

    // === Create capsules ===
    let mut capsules: Vec<(BodyHandle, f32, f32, Color)> = Vec::new();

    let capsule_configs = [
        (PhysVec3::new(-6.0, 8.0, 0.0), 0.2, 0.5, PINK),
        (PhysVec3::new(6.0, 8.0, 0.0), 0.2, 0.5, PINK),
        (PhysVec3::new(0.0, 20.0, 0.0), 0.25, 0.6, LIME),
    ];

    for (pos, radius, half_height, color) in capsule_configs {
        let capsule = world.create_body(
            RigidBodyDesc::dynamic()
                .with_position(pos)
                .with_mass(1.5),
        );
        world.attach_collider(capsule, Shape::capsule(radius, half_height), 1.0);
        capsules.push((capsule, radius, half_height, color));
    }

    // Track mouse-spawned bodies with their shape info
    enum SpawnedShape {
        Ball(f32, Color),
        Box(PhysVec3, Color),
        Capsule(f32, f32, Color),
    }
    let mut spawned_bodies: Vec<(BodyHandle, SpawnedShape)> = Vec::new();
    let mut spawn_mode = 0; // 0 = ball, 1 = box, 2 = capsule

    let dt = 1.0 / 60.0;
    let mut paused = false;
    let mut time_scale = 1.0; // Normal speed by default

    loop {
        // Handle input
        if is_key_pressed(KeyCode::Space) {
            paused = !paused;
        }

        if is_key_pressed(KeyCode::Key1) {
            spawn_mode = 0; // Ball
        }
        if is_key_pressed(KeyCode::Key2) {
            spawn_mode = 1; // Box
        }
        if is_key_pressed(KeyCode::Key3) {
            spawn_mode = 2; // Capsule
        }

        // Time scale controls
        if is_key_pressed(KeyCode::Up) {
            time_scale = (time_scale * 1.5_f32).min(2.0);
        }
        if is_key_pressed(KeyCode::Down) {
            time_scale = (time_scale / 1.5_f32).max(0.05);
        }

        if is_key_pressed(KeyCode::R) {
            // Reset - remove spawned bodies
            for (handle, _) in spawned_bodies.drain(..) {
                world.remove_body(handle);
            }
            // Reset dynamic objects to initial positions
            for (i, (ball, _, _)) in balls.iter().enumerate() {
                if let Some(body) = world.body_mut(*ball) {
                    if i < 3 {
                        body.position = PhysVec3::new((i as f32 - 1.0) * 2.0, 12.0 + i as f32 * 1.0, 0.0);
                    } else if i < 7 {
                        let j = i - 3;
                        body.position = PhysVec3::new((j as f32 - 1.5) * 1.5, 15.0 + j as f32 * 0.8, 0.0);
                    } else {
                        let j = i - 7;
                        body.position = PhysVec3::new((j as f32 - 2.5) * 1.0, 18.0 + j as f32 * 0.5, 0.0);
                    }
                    body.linear_velocity = PhysVec3::ZERO;
                    body.angular_velocity = PhysVec3::ZERO;
                    body.rotation = PhysQuat::IDENTITY;
                    body.wake_up();
                }
            }
            for (i, (box_body, _, _)) in boxes.iter().enumerate() {
                if let Some(body) = world.body_mut(*box_body) {
                    let positions = [
                        PhysVec3::new(-4.0, 14.0, 0.0),
                        PhysVec3::new(4.0, 14.0, 0.0),
                        PhysVec3::new(0.0, 16.0, 0.0),
                        PhysVec3::new(-2.0, 17.0, 0.0),
                        PhysVec3::new(2.0, 17.0, 0.0),
                    ];
                    if i < positions.len() {
                        body.position = positions[i];
                    }
                    body.linear_velocity = PhysVec3::ZERO;
                    body.angular_velocity = PhysVec3::ZERO;
                    body.rotation = PhysQuat::IDENTITY;
                    body.wake_up();
                }
            }
            for (i, (capsule, _, _, _)) in capsules.iter().enumerate() {
                if let Some(body) = world.body_mut(*capsule) {
                    let positions = [
                        PhysVec3::new(-6.0, 8.0, 0.0),
                        PhysVec3::new(6.0, 8.0, 0.0),
                        PhysVec3::new(0.0, 20.0, 0.0),
                    ];
                    if i < positions.len() {
                        body.position = positions[i];
                    }
                    body.linear_velocity = PhysVec3::ZERO;
                    body.angular_velocity = PhysVec3::ZERO;
                    body.rotation = PhysQuat::IDENTITY;
                    body.wake_up();
                }
            }
        }

        // Spawn on mouse click
        if is_mouse_button_pressed(MouseButton::Left) {
            let (mx, my) = mouse_position();
            let px = (mx - WINDOW_WIDTH / 2.0) / SCALE;
            let py = (WINDOW_HEIGHT - 80.0 - my) / SCALE;

            match spawn_mode {
                0 => {
                    // Ball
                    let radius = 0.3;
                    let new_ball = world.create_body(
                        RigidBodyDesc::dynamic()
                            .with_position(PhysVec3::new(px, py, 0.0))
                            .with_mass(1.0),
                    );
                    world.attach_collider(new_ball, Shape::sphere(radius), 1.0);
                    spawned_bodies.push((new_ball, SpawnedShape::Ball(radius, GOLD)));
                }
                1 => {
                    // Box
                    let extents = PhysVec3::new(0.35, 0.35, 0.35);
                    let new_box = world.create_body(
                        RigidBodyDesc::dynamic()
                            .with_position(PhysVec3::new(px, py, 0.0))
                            .with_mass(1.5),
                    );
                    world.attach_collider(new_box, Shape::cuboid(extents), 1.0);
                    spawned_bodies.push((new_box, SpawnedShape::Box(extents, GREEN)));
                }
                2 => {
                    // Capsule
                    let radius = 0.2;
                    let half_height = 0.4;
                    let new_capsule = world.create_body(
                        RigidBodyDesc::dynamic()
                            .with_position(PhysVec3::new(px, py, 0.0))
                            .with_mass(1.2),
                    );
                    world.attach_collider(new_capsule, Shape::capsule(radius, half_height), 1.0);
                    spawned_bodies.push((new_capsule, SpawnedShape::Capsule(radius, half_height, MAGENTA)));
                }
                _ => {}
            }
        }

        // Step physics with time scale
        if !paused {
            world.step(dt * time_scale);
        }

        // === DRAWING ===
        clear_background(Color::from_rgba(25, 25, 35, 255));

        // Draw platforms
        for platform in &platforms {
            let pos = world.body_position(platform.handle);
            let rot = world.body_rotation(platform.handle);
            let (sx, sy) = physics_to_screen(pos);
            let (_, _, angle) = rot.to_euler();

            draw_rectangle_ex(
                sx,
                sy,
                platform.half_extents.x * 2.0 * SCALE,
                platform.half_extents.y * 2.0 * SCALE,
                DrawRectangleParams {
                    offset: macroquad::math::Vec2::new(0.5, 0.5),
                    rotation: -angle,
                    color: platform.color,
                },
            );
        }

        // Draw balls
        for (ball, radius, color) in &balls {
            let pos = world.body_position(*ball);
            let (sx, sy) = physics_to_screen(pos);
            draw_circle(sx, sy, radius * SCALE, *color);
            draw_circle_lines(sx, sy, radius * SCALE, 2.0, Color::from_rgba(50, 50, 50, 255));
        }

        // Draw boxes
        for (box_body, extents, color) in &boxes {
            let pos = world.body_position(*box_body);
            let rot = world.body_rotation(*box_body);
            let (sx, sy) = physics_to_screen(pos);
            let (_, _, angle) = rot.to_euler();

            draw_rectangle_ex(
                sx,
                sy,
                extents.x * 2.0 * SCALE,
                extents.y * 2.0 * SCALE,
                DrawRectangleParams {
                    offset: macroquad::math::Vec2::new(0.5, 0.5),
                    rotation: -angle,
                    color: *color,
                },
            );
        }

        // Draw capsules (as rounded rectangles approximation)
        for (capsule, radius, half_height, color) in &capsules {
            let pos = world.body_position(*capsule);
            let rot = world.body_rotation(*capsule);
            let (sx, sy) = physics_to_screen(pos);
            let (_, _, angle) = rot.to_euler();

            // Draw capsule as two circles + rectangle
            let r = radius * SCALE;
            let h = half_height * SCALE;

            // Calculate rotated offset for the two end circles
            let cos_a = angle.cos();
            let sin_a = angle.sin();
            let dy_x = -h * sin_a;
            let dy_y = -h * cos_a;

            // Top circle
            draw_circle(sx + dy_x, sy + dy_y, r, *color);
            // Bottom circle
            draw_circle(sx - dy_x, sy - dy_y, r, *color);
            // Middle rectangle
            draw_rectangle_ex(
                sx,
                sy,
                r * 2.0,
                h * 2.0,
                DrawRectangleParams {
                    offset: macroquad::math::Vec2::new(0.5, 0.5),
                    rotation: -angle,
                    color: *color,
                },
            );
        }

        // Draw spawned bodies
        for (handle, shape) in &spawned_bodies {
            if let Some(body) = world.body(*handle) {
                let pos = body.position;
                let rot = body.rotation;
                let (sx, sy) = physics_to_screen(pos);
                let (_, _, angle) = rot.to_euler();

                match shape {
                    SpawnedShape::Ball(radius, color) => {
                        draw_circle(sx, sy, radius * SCALE, *color);
                        draw_circle_lines(sx, sy, radius * SCALE, 2.0, WHITE);
                    }
                    SpawnedShape::Box(extents, color) => {
                        draw_rectangle_ex(
                            sx,
                            sy,
                            extents.x * 2.0 * SCALE,
                            extents.y * 2.0 * SCALE,
                            DrawRectangleParams {
                                offset: macroquad::math::Vec2::new(0.5, 0.5),
                                rotation: -angle,
                                color: *color,
                            },
                        );
                    }
                    SpawnedShape::Capsule(radius, half_height, color) => {
                        let r = radius * SCALE;
                        let h = half_height * SCALE;
                        let cos_a = angle.cos();
                        let sin_a = angle.sin();
                        let dy_x = -h * sin_a;
                        let dy_y = -h * cos_a;

                        draw_circle(sx + dy_x, sy + dy_y, r, *color);
                        draw_circle(sx - dy_x, sy - dy_y, r, *color);
                        draw_rectangle_ex(
                            sx,
                            sy,
                            r * 2.0,
                            h * 2.0,
                            DrawRectangleParams {
                                offset: macroquad::math::Vec2::new(0.5, 0.5),
                                rotation: -angle,
                                color: *color,
                            },
                        );
                    }
                }
            }
        }

        // === UI ===
        draw_text("RustPhy Visual Demo", 10.0, 25.0, 28.0, WHITE);
        draw_text(
            &format!("Bodies: {}", world.num_bodies()),
            10.0,
            50.0,
            20.0,
            LIGHTGRAY,
        );
        draw_text(
            &format!("FPS: {}", get_fps()),
            10.0,
            70.0,
            20.0,
            LIGHTGRAY,
        );
        draw_text(
            &format!("Speed: {:.0}%", time_scale * 100.0),
            10.0,
            90.0,
            20.0,
            LIGHTGRAY,
        );

        // Spawn mode indicator
        let mode_text = match spawn_mode {
            0 => "Spawn: [1] BALL",
            1 => "Spawn: [2] BOX",
            2 => "Spawn: [3] CAPSULE",
            _ => "",
        };
        draw_text(mode_text, 10.0, 115.0, 20.0, GOLD);

        if paused {
            draw_text("PAUSED", WINDOW_WIDTH / 2.0 - 60.0, 35.0, 36.0, YELLOW);
        }

        draw_text(
            "Controls: [Space] Pause | [R] Reset | [1/2/3] Shape | [Up/Down] Speed | [Click] Spawn",
            10.0,
            WINDOW_HEIGHT - 10.0,
            16.0,
            GRAY,
        );

        next_frame().await
    }
}
