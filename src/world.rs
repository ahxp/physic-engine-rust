use std::collections::HashMap;

use crate::collision::{
    epa, gjk, BodyHandle, Bvh, CollisionPair, ContactManifold, ContactPoint, GjkResult,
};
use crate::dynamics::{integrate_positions, integrate_velocities, BodyType, RigidBody, RigidBodyDesc};
use crate::geometry::{MassProperties, Shape};
use crate::math::{Quat, Transform, Vec3};
use crate::solver::{solve_position_constraints, PgsSolver, SolverConfig};

/// Configuration for the physics world
#[derive(Debug, Clone)]
pub struct WorldConfig {
    /// Gravity vector
    pub gravity: Vec3,
    /// Solver configuration
    pub solver: SolverConfig,
    /// Maximum substeps per frame
    pub max_substeps: usize,
    /// Sleep threshold (kinetic energy)
    pub sleep_threshold: f32,
    /// Time before a body can sleep
    pub sleep_time_threshold: f32,
}

impl Default for WorldConfig {
    fn default() -> Self {
        Self {
            gravity: Vec3::new(0.0, -9.81, 0.0),
            solver: SolverConfig::default(),
            max_substeps: 8,           // More substeps for stability
            sleep_threshold: 0.1,      // Kinetic energy threshold for sleep
            sleep_time_threshold: 0.3, // Faster sleep for stability
        }
    }
}

/// The main physics world containing all bodies and managing simulation
pub struct World {
    /// Configuration
    config: WorldConfig,
    /// All rigid bodies
    bodies: Vec<RigidBody>,
    /// Free body indices for reuse
    free_bodies: Vec<usize>,
    /// Broad phase collision detection
    broad_phase: Bvh,
    /// Contact manifolds indexed by collision pair
    manifolds: HashMap<CollisionPair, ContactManifold>,
    /// Active manifolds for current step
    active_manifolds: Vec<ContactManifold>,
    /// Constraint solver
    solver: PgsSolver,
    /// Current simulation time
    time: f32,
}

impl Default for World {
    fn default() -> Self {
        Self::new(WorldConfig::default())
    }
}

impl World {
    /// Creates a new physics world with the given configuration
    pub fn new(config: WorldConfig) -> Self {
        Self {
            solver: PgsSolver::new(config.solver.clone()),
            config,
            bodies: Vec::new(),
            free_bodies: Vec::new(),
            broad_phase: Bvh::new(),
            manifolds: HashMap::new(),
            active_manifolds: Vec::new(),
            time: 0.0,
        }
    }

    /// Creates a new rigid body and returns its handle
    pub fn create_body(&mut self, desc: RigidBodyDesc) -> BodyHandle {
        let handle = if let Some(index) = self.free_bodies.pop() {
            BodyHandle::new(index as u32)
        } else {
            let index = self.bodies.len();
            self.bodies.push(RigidBody::default());
            BodyHandle::new(index as u32)
        };

        let body = &mut self.bodies[handle.index()];
        body.handle = handle;
        body.body_type = desc.body_type;
        body.position = desc.position;
        body.rotation = desc.rotation;
        body.linear_velocity = desc.linear_velocity;
        body.angular_velocity = desc.angular_velocity;
        body.friction = desc.friction;
        body.restitution = desc.restitution;
        body.linear_damping = desc.linear_damping;
        body.angular_damping = desc.angular_damping;
        body.can_sleep = desc.can_sleep;
        body.user_data = desc.user_data;

        if desc.body_type == BodyType::Dynamic && desc.mass > 0.0 {
            body.inv_mass = 1.0 / desc.mass;
        } else {
            body.inv_mass = 0.0;
        }

        body.update_world_inertia();

        handle
    }

    /// Attaches a collision shape to a body
    pub fn attach_collider(&mut self, handle: BodyHandle, shape: Shape, density: f32) {
        let index = handle.index();
        if index >= self.bodies.len() {
            return;
        }

        let body = &mut self.bodies[index];
        body.shape = Some(shape);

        // Update mass properties from shape
        if body.is_dynamic() {
            let props = shape.mass_properties(density);
            body.inv_mass = props.inv_mass();
            body.inv_inertia_local = props.inv_inertia();
            body.update_world_inertia();
        }

        // Add to broad phase
        let aabb = shape.world_aabb(body.transform());
        self.broad_phase.insert(handle, aabb);
    }

    /// Removes a body from the world
    pub fn remove_body(&mut self, handle: BodyHandle) {
        let index = handle.index();
        if index >= self.bodies.len() {
            return;
        }

        self.broad_phase.remove(handle);
        self.bodies[index] = RigidBody::default();
        self.free_bodies.push(index);

        // Remove any manifolds involving this body
        self.manifolds.retain(|pair, _| pair.body_a != handle && pair.body_b != handle);
    }

    /// Gets a reference to a body
    pub fn body(&self, handle: BodyHandle) -> Option<&RigidBody> {
        self.bodies.get(handle.index())
    }

    /// Gets a mutable reference to a body
    pub fn body_mut(&mut self, handle: BodyHandle) -> Option<&mut RigidBody> {
        self.bodies.get_mut(handle.index())
    }

    /// Gets the position of a body
    pub fn body_position(&self, handle: BodyHandle) -> Vec3 {
        self.bodies
            .get(handle.index())
            .map(|b| b.position)
            .unwrap_or(Vec3::ZERO)
    }

    /// Gets the rotation of a body
    pub fn body_rotation(&self, handle: BodyHandle) -> Quat {
        self.bodies
            .get(handle.index())
            .map(|b| b.rotation)
            .unwrap_or(Quat::IDENTITY)
    }

    /// Gets the transform of a body
    pub fn body_transform(&self, handle: BodyHandle) -> Transform {
        self.bodies
            .get(handle.index())
            .map(|b| b.transform())
            .unwrap_or(Transform::IDENTITY)
    }

    /// Sets the position of a body
    pub fn set_body_position(&mut self, handle: BodyHandle, position: Vec3) {
        if let Some(body) = self.bodies.get_mut(handle.index()) {
            body.position = position;
            body.wake_up();
            self.update_body_aabb(handle);
        }
    }

    /// Sets the rotation of a body
    pub fn set_body_rotation(&mut self, handle: BodyHandle, rotation: Quat) {
        if let Some(body) = self.bodies.get_mut(handle.index()) {
            body.rotation = rotation;
            body.update_world_inertia();
            body.wake_up();
            self.update_body_aabb(handle);
        }
    }

    /// Sets the linear velocity of a body
    pub fn set_linear_velocity(&mut self, handle: BodyHandle, velocity: Vec3) {
        if let Some(body) = self.bodies.get_mut(handle.index()) {
            body.linear_velocity = velocity;
            body.wake_up();
        }
    }

    /// Sets the angular velocity of a body
    pub fn set_angular_velocity(&mut self, handle: BodyHandle, velocity: Vec3) {
        if let Some(body) = self.bodies.get_mut(handle.index()) {
            body.angular_velocity = velocity;
            body.wake_up();
        }
    }

    /// Applies a force to a body at its center of mass
    pub fn apply_force(&mut self, handle: BodyHandle, force: Vec3) {
        if let Some(body) = self.bodies.get_mut(handle.index()) {
            body.apply_force(force);
        }
    }

    /// Applies a force to a body at a world point
    pub fn apply_force_at_point(&mut self, handle: BodyHandle, force: Vec3, point: Vec3) {
        if let Some(body) = self.bodies.get_mut(handle.index()) {
            body.apply_force_at_point(force, point);
        }
    }

    /// Applies an impulse to a body at its center of mass
    pub fn apply_impulse(&mut self, handle: BodyHandle, impulse: Vec3) {
        if let Some(body) = self.bodies.get_mut(handle.index()) {
            body.apply_impulse(impulse);
        }
    }

    /// Applies an impulse to a body at a world point
    pub fn apply_impulse_at_point(&mut self, handle: BodyHandle, impulse: Vec3, point: Vec3) {
        if let Some(body) = self.bodies.get_mut(handle.index()) {
            body.apply_impulse_at_point(impulse, point);
        }
    }

    /// Sets the gravity
    pub fn set_gravity(&mut self, gravity: Vec3) {
        self.config.gravity = gravity;
    }

    /// Gets the gravity
    pub fn gravity(&self) -> Vec3 {
        self.config.gravity
    }

    /// Steps the simulation by the given time delta
    pub fn step(&mut self, dt: f32) {
        if dt <= 0.0 {
            return;
        }

        // Use fixed timestep with sub-stepping for stability
        const FIXED_DT: f32 = 1.0 / 60.0; // 60 Hz internal simulation
        let max_substeps = self.config.max_substeps.max(1);

        let num_substeps = ((dt / FIXED_DT).ceil() as usize).min(max_substeps);
        let substep_dt = dt / num_substeps as f32;

        for _ in 0..num_substeps {
            self.substep(substep_dt);
        }

        self.time += dt;
    }

    /// Performs a single simulation substep
    fn substep(&mut self, dt: f32) {
        // Integrate velocities (apply gravity and forces)
        self.integrate_velocities(dt);

        // Detect collisions
        self.detect_collisions();

        // Prepare constraints
        self.solver.prepare(&self.active_manifolds, &self.bodies, dt);

        // Warm start
        self.solver.warm_start(&mut self.bodies);

        // Solve velocity constraints
        self.solver.solve_velocity(&mut self.bodies);

        // Integrate positions
        self.integrate_positions(dt);

        // Solve position constraints to resolve remaining penetration
        solve_position_constraints(&self.active_manifolds, &mut self.bodies, self.solver.config());

        // Store impulses for warm starting next frame
        self.solver.store_impulses(&mut self.active_manifolds);

        // Update AABBs and sleep state
        self.post_step(dt);
    }

    /// Integrates velocities for all bodies
    fn integrate_velocities(&mut self, dt: f32) {
        for body in &mut self.bodies {
            if body.handle.is_valid() {
                integrate_velocities(body, self.config.gravity, dt);
            }
        }
    }

    /// Integrates positions for all bodies
    fn integrate_positions(&mut self, dt: f32) {
        for body in &mut self.bodies {
            if body.handle.is_valid() {
                integrate_positions(body, dt);
            }
        }
    }

    /// Detects collisions and generates contact manifolds
    fn detect_collisions(&mut self) {
        self.active_manifolds.clear();

        // Broad phase
        let pairs = self.broad_phase.query_pairs();

        // Narrow phase
        for (handle_a, handle_b) in pairs {
            let body_a = &self.bodies[handle_a.index()];
            let body_b = &self.bodies[handle_b.index()];

            // Skip if both are static or sleeping
            if (!body_a.is_dynamic() && !body_b.is_dynamic())
                || (!body_a.is_awake && !body_b.is_awake)
            {
                continue;
            }

            let shape_a = match &body_a.shape {
                Some(s) => s,
                None => continue,
            };
            let shape_b = match &body_b.shape {
                Some(s) => s,
                None => continue,
            };

            let transform_a = body_a.transform();
            let transform_b = body_b.transform();

            // GJK collision detection
            let result = gjk(shape_a, transform_a, shape_b, transform_b);

            if let GjkResult::Intersecting(simplex) = result {
                // EPA to find penetration
                if let Some(epa_result) = epa(&simplex, shape_a, transform_a, shape_b, transform_b) {
                    let pair = CollisionPair::new(handle_a, handle_b);

                    // Get or create manifold
                    let old_manifold = self.manifolds.get(&pair).cloned();
                    let manifold = self.manifolds.entry(pair).or_insert_with(|| {
                        let mut m = ContactManifold::new(handle_a, handle_b);
                        m.friction = (body_a.friction * body_b.friction).sqrt();
                        m.restitution = (body_a.restitution + body_b.restitution) * 0.5;
                        m
                    });

                    // Clear old contacts and create new one
                    manifold.clear();

                    let contact = ContactPoint::new(
                        epa_result.point_a,
                        epa_result.point_b,
                        transform_a.inverse_transform_point(epa_result.point_a),
                        transform_b.inverse_transform_point(epa_result.point_b),
                        epa_result.normal,
                        epa_result.depth,
                    );

                    manifold.add_point(contact);

                    // Warm start from previous manifold
                    if let Some(old) = old_manifold {
                        manifold.warm_start(&old);
                    }

                    self.active_manifolds.push(manifold.clone());

                    // Wake up bodies
                    self.bodies[handle_a.index()].wake_up();
                    self.bodies[handle_b.index()].wake_up();
                }
            }
        }

        // Clean up stale manifolds
        let active_pairs: std::collections::HashSet<_> = self
            .active_manifolds
            .iter()
            .map(|m| CollisionPair::new(m.body_a, m.body_b))
            .collect();
        self.manifolds.retain(|pair, _| active_pairs.contains(pair));
    }

    /// Post-step updates (AABBs, sleep)
    fn post_step(&mut self, dt: f32) {
        for body in &mut self.bodies {
            if !body.handle.is_valid() {
                continue;
            }

            // Update AABB
            if let Some(shape) = &body.shape {
                let aabb = shape.world_aabb(body.transform());
                self.broad_phase.update(body.handle, aabb);
            }

            // Update sleep state
            body.update_sleep(dt, self.config.sleep_threshold, self.config.sleep_time_threshold);

            // Clear forces
            body.clear_forces();
        }
    }

    /// Updates the AABB of a body in the broad phase
    fn update_body_aabb(&mut self, handle: BodyHandle) {
        if let Some(body) = self.bodies.get(handle.index()) {
            if let Some(shape) = &body.shape {
                let aabb = shape.world_aabb(body.transform());
                self.broad_phase.update(handle, aabb);
            }
        }
    }

    /// Returns the number of bodies in the world
    pub fn num_bodies(&self) -> usize {
        self.bodies.len() - self.free_bodies.len()
    }

    /// Returns the current simulation time
    pub fn time(&self) -> f32 {
        self.time
    }

    /// Returns an iterator over all body handles
    pub fn bodies(&self) -> impl Iterator<Item = BodyHandle> + '_ {
        self.bodies
            .iter()
            .filter(|b| b.handle.is_valid())
            .map(|b| b.handle)
    }

    /// Performs a ray cast and returns the first hit
    pub fn ray_cast(&self, origin: Vec3, direction: Vec3, max_distance: f32) -> Option<RayCastHit> {
        let candidates = self.broad_phase.query_ray(origin, direction, max_distance);

        let mut closest: Option<RayCastHit> = None;

        for handle in candidates {
            let body = &self.bodies[handle.index()];
            let shape = match &body.shape {
                Some(s) => s,
                None => continue,
            };

            // Simple sphere/box intersection for now
            // A proper implementation would use shape-specific ray casting
            let aabb = shape.world_aabb(body.transform());
            if let Some((t, _)) = aabb.ray_intersection(origin, direction) {
                if t <= max_distance {
                    let hit = RayCastHit {
                        body: handle,
                        point: origin + direction * t,
                        normal: Vec3::Y, // Simplified
                        distance: t,
                    };

                    if closest.is_none() || t < closest.as_ref().unwrap().distance {
                        closest = Some(hit);
                    }
                }
            }
        }

        closest
    }
}

/// Result of a ray cast query
#[derive(Debug, Clone, Copy)]
pub struct RayCastHit {
    /// Body that was hit
    pub body: BodyHandle,
    /// World space hit point
    pub point: Vec3,
    /// Surface normal at hit point
    pub normal: Vec3,
    /// Distance from ray origin
    pub distance: f32,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Shape;

    #[test]
    fn test_world_creation() {
        let world = World::default();
        assert_eq!(world.num_bodies(), 0);
    }

    #[test]
    fn test_create_body() {
        let mut world = World::default();

        let handle = world.create_body(RigidBodyDesc::dynamic().with_position(Vec3::new(0.0, 5.0, 0.0)));

        assert_eq!(world.num_bodies(), 1);
        assert_eq!(world.body_position(handle), Vec3::new(0.0, 5.0, 0.0));
    }

    #[test]
    fn test_gravity_simulation() {
        let mut world = World::default();
        world.set_gravity(Vec3::new(0.0, -10.0, 0.0));

        let handle = world.create_body(
            RigidBodyDesc::dynamic()
                .with_position(Vec3::new(0.0, 10.0, 0.0))
                .with_mass(1.0),
        );
        world.attach_collider(handle, Shape::sphere(1.0), 1.0);

        // Step simulation
        for _ in 0..60 {
            world.step(1.0 / 60.0);
        }

        // Body should have fallen
        let pos = world.body_position(handle);
        assert!(pos.y < 10.0);
    }

    #[test]
    fn test_collision_detection() {
        let mut world = World::default();
        world.set_gravity(Vec3::new(0.0, -10.0, 0.0));

        // Create floor at Y=0 with half-height 0.5 (top surface at Y=0.5)
        let floor = world.create_body(RigidBodyDesc::fixed().with_position(Vec3::ZERO));
        world.attach_collider(floor, Shape::cuboid(Vec3::new(10.0, 0.5, 10.0)), 1.0);

        // Create ball above floor at Y=3
        let ball = world.create_body(
            RigidBodyDesc::dynamic()
                .with_position(Vec3::new(0.0, 3.0, 0.0))
                .with_mass(1.0),
        );
        world.attach_collider(ball, Shape::sphere(0.5), 1.0);

        // Step simulation
        for _ in 0..300 {
            world.step(1.0 / 60.0);
        }

        // Ball should have fallen and be resting on floor
        // Floor top is at Y=0.5, ball radius is 0.5, so ball center should be around Y=1.0
        let pos = world.body_position(ball);
        assert!(pos.y > 0.5, "Ball fell through floor: y={}", pos.y);
        assert!(pos.y < 3.0, "Ball didn't fall: y={}", pos.y);
    }

    #[test]
    fn test_remove_body() {
        let mut world = World::default();

        let handle = world.create_body(RigidBodyDesc::dynamic());
        assert_eq!(world.num_bodies(), 1);

        world.remove_body(handle);
        assert_eq!(world.num_bodies(), 0);
    }

    #[test]
    fn test_box_stack() {
        let mut world = World::default();
        world.set_gravity(Vec3::new(0.0, -10.0, 0.0));

        // Create floor
        let floor = world.create_body(RigidBodyDesc::fixed().with_position(Vec3::ZERO));
        world.attach_collider(floor, Shape::cuboid(Vec3::new(10.0, 0.5, 10.0)), 1.0);

        // Create stack of 3 boxes
        let box_size = 0.5; // half-extent
        let mut boxes = Vec::new();
        for i in 0..3 {
            let y = 0.5 + box_size + (i as f32) * (box_size * 2.0 + 0.1); // Small gap
            let handle = world.create_body(
                RigidBodyDesc::dynamic()
                    .with_position(Vec3::new(0.0, y, 0.0))
                    .with_mass(1.0),
            );
            world.attach_collider(handle, Shape::cuboid(Vec3::new(box_size, box_size, box_size)), 1.0);
            boxes.push(handle);
        }

        // Step simulation until settled
        for _ in 0..600 {
            world.step(1.0 / 60.0);
        }

        // All boxes should be stacked on floor
        // Floor top at 0.5, box half-height 0.5
        // Box 0 center should be at ~1.0
        // Box 1 center should be at ~2.0
        // Box 2 center should be at ~3.0
        let pos0 = world.body_position(boxes[0]);
        let pos1 = world.body_position(boxes[1]);
        let pos2 = world.body_position(boxes[2]);

        assert!(pos0.y > 0.8 && pos0.y < 1.5, "Box 0 position wrong: y={}", pos0.y);
        assert!(pos1.y > 1.8 && pos1.y < 2.5, "Box 1 position wrong: y={}", pos1.y);
        assert!(pos2.y > 2.8 && pos2.y < 3.5, "Box 2 position wrong: y={}", pos2.y);

        // Check they're roughly stacked (small x deviation)
        assert!(pos0.x.abs() < 0.5, "Box 0 drifted: x={}", pos0.x);
        assert!(pos1.x.abs() < 0.5, "Box 1 drifted: x={}", pos1.x);
        assert!(pos2.x.abs() < 0.5, "Box 2 drifted: x={}", pos2.x);
    }
}
