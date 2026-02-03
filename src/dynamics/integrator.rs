use crate::math::Vec3;

use super::rigid_body::RigidBody;

/// Integration methods for physics simulation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntegrationMethod {
    /// Semi-implicit Euler (stable, fast)
    SemiImplicitEuler,
    /// Explicit Euler (simple, less stable)
    ExplicitEuler,
}

impl Default for IntegrationMethod {
    fn default() -> Self {
        IntegrationMethod::SemiImplicitEuler
    }
}

/// Integrates velocities (applies forces to velocities)
pub fn integrate_velocities(body: &mut RigidBody, gravity: Vec3, dt: f32) {
    if !body.is_dynamic() || !body.is_awake {
        return;
    }

    // Apply gravity
    body.linear_velocity += gravity * dt;

    // Apply accumulated forces
    body.linear_velocity += body.force * body.inv_mass * dt;
    body.angular_velocity += body.inv_inertia_world * body.torque * dt;

    // Apply damping
    body.linear_velocity *= (1.0 - body.linear_damping).powf(dt);
    body.angular_velocity *= (1.0 - body.angular_damping).powf(dt);

    // Clamp velocities to prevent instability
    const MAX_LINEAR_VELOCITY: f32 = 100.0;
    const MAX_ANGULAR_VELOCITY: f32 = 50.0;

    let linear_speed = body.linear_velocity.length();
    if linear_speed > MAX_LINEAR_VELOCITY {
        body.linear_velocity *= MAX_LINEAR_VELOCITY / linear_speed;
    }

    let angular_speed = body.angular_velocity.length();
    if angular_speed > MAX_ANGULAR_VELOCITY {
        body.angular_velocity *= MAX_ANGULAR_VELOCITY / angular_speed;
    }
}

/// Integrates positions (applies velocities to positions)
pub fn integrate_positions(body: &mut RigidBody, dt: f32) {
    if !body.is_dynamic() || !body.is_awake {
        return;
    }

    // Update position
    body.position += body.linear_velocity * dt;

    // Update rotation using quaternion integration
    body.rotation = body.rotation.integrate(body.angular_velocity, dt);

    // Update world inertia tensor
    body.update_world_inertia();
}

/// Performs a full integration step (semi-implicit Euler)
pub fn integrate_semi_implicit_euler(body: &mut RigidBody, gravity: Vec3, dt: f32) {
    // First update velocities (using current forces)
    integrate_velocities(body, gravity, dt);

    // Then update positions (using new velocities)
    integrate_positions(body, dt);

    // Clear forces for next frame
    body.clear_forces();
}

/// Performs a full integration step (explicit Euler)
pub fn integrate_explicit_euler(body: &mut RigidBody, gravity: Vec3, dt: f32) {
    if !body.is_dynamic() || !body.is_awake {
        return;
    }

    // Save current velocities
    let v = body.linear_velocity;
    let w = body.angular_velocity;

    // Update velocities
    integrate_velocities(body, gravity, dt);

    // Update positions using old velocities
    body.position += v * dt;
    body.rotation = body.rotation.integrate(w, dt);

    // Update world inertia tensor
    body.update_world_inertia();

    // Clear forces
    body.clear_forces();
}

/// Solves position constraints by directly moving bodies
pub fn solve_position_constraint(
    body_a: &mut RigidBody,
    body_b: &mut RigidBody,
    contact_point: Vec3,
    normal: Vec3,
    penetration: f32,
    baumgarte: f32,
) {
    if penetration <= 0.0 {
        return;
    }

    let correction = penetration * baumgarte;

    let r_a = contact_point - body_a.position;
    let r_b = contact_point - body_b.position;

    // Compute effective mass
    let inv_mass_a = body_a.inv_mass;
    let inv_mass_b = body_b.inv_mass;

    let rn_a = r_a.cross(normal);
    let rn_b = r_b.cross(normal);

    let k = inv_mass_a
        + inv_mass_b
        + rn_a.dot(body_a.inv_inertia_world * rn_a)
        + rn_b.dot(body_b.inv_inertia_world * rn_b);

    if k <= 0.0 {
        return;
    }

    let impulse = correction / k;

    // Apply position correction
    body_a.position += normal * (impulse * inv_mass_a);
    body_b.position -= normal * (impulse * inv_mass_b);

    // Apply rotation correction
    if body_a.is_dynamic() {
        let angular_impulse = body_a.inv_inertia_world * rn_a * impulse;
        body_a.rotation = body_a.rotation.integrate(angular_impulse, 1.0);
        body_a.update_world_inertia();
    }

    if body_b.is_dynamic() {
        let angular_impulse = body_b.inv_inertia_world * rn_b * impulse;
        body_b.rotation = body_b.rotation.integrate(-angular_impulse, 1.0);
        body_b.update_world_inertia();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Quat;

    #[test]
    fn test_gravity_integration() {
        let mut body = RigidBody::new().with_mass(1.0).with_position(Vec3::ZERO);

        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let dt = 1.0 / 60.0;

        integrate_semi_implicit_euler(&mut body, gravity, dt);

        // Velocity should have increased due to gravity
        assert!(body.linear_velocity.y < 0.0);

        // Position should have changed
        assert!(body.position.y < 0.0);
    }

    #[test]
    fn test_static_body_no_integration() {
        use super::super::rigid_body::BodyType;

        let mut body = RigidBody::new()
            .with_type(BodyType::Static)
            .with_position(Vec3::ZERO);

        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let dt = 1.0 / 60.0;

        integrate_semi_implicit_euler(&mut body, gravity, dt);

        // Position should not change for static body
        assert_eq!(body.position, Vec3::ZERO);
        assert_eq!(body.linear_velocity, Vec3::ZERO);
    }

    #[test]
    fn test_angular_velocity_integration() {
        let mut body = RigidBody::new()
            .with_mass(1.0)
            .with_position(Vec3::ZERO)
            .with_rotation(Quat::IDENTITY);

        body.angular_velocity = Vec3::new(0.0, 0.0, std::f32::consts::PI); // 180 deg/s around Z

        let gravity = Vec3::ZERO;
        let dt = 1.0; // 1 second

        integrate_semi_implicit_euler(&mut body, gravity, dt);

        // After 1 second at PI rad/s, should have rotated ~180 degrees
        // The local X axis should now point in -X direction
        let local_x = body.rotation.local_x();
        assert!(local_x.x < -0.9);
    }

    #[test]
    fn test_damping() {
        let mut body = RigidBody::new()
            .with_mass(1.0)
            .with_position(Vec3::ZERO)
            .with_linear_damping(0.1);

        body.linear_velocity = Vec3::new(10.0, 0.0, 0.0);

        let gravity = Vec3::ZERO;
        let dt = 1.0;

        integrate_semi_implicit_euler(&mut body, gravity, dt);

        // Velocity should have decreased due to damping
        assert!(body.linear_velocity.x < 10.0);
        assert!(body.linear_velocity.x > 0.0);
    }
}
