use crate::collision::ContactPoint;
use crate::dynamics::RigidBody;
use crate::math::Vec3;

/// A velocity constraint for a contact point
#[derive(Debug, Clone, Copy)]
pub struct ContactVelocityConstraint {
    /// Contact point in world space
    pub point: Vec3,
    /// Contact normal (pointing from B to A)
    pub normal: Vec3,
    /// First tangent direction
    pub tangent1: Vec3,
    /// Second tangent direction
    pub tangent2: Vec3,
    /// Radius vector from body A center to contact point
    pub r_a: Vec3,
    /// Radius vector from body B center to contact point
    pub r_b: Vec3,
    /// Effective mass for normal constraint
    pub normal_mass: f32,
    /// Effective mass for tangent1 constraint
    pub tangent1_mass: f32,
    /// Effective mass for tangent2 constraint
    pub tangent2_mass: f32,
    /// Velocity bias for restitution
    pub velocity_bias: f32,
    /// Accumulated normal impulse
    pub normal_impulse: f32,
    /// Accumulated tangent1 impulse
    pub tangent1_impulse: f32,
    /// Accumulated tangent2 impulse
    pub tangent2_impulse: f32,
}

impl ContactVelocityConstraint {
    /// Creates a new contact velocity constraint
    pub fn new(
        contact: &ContactPoint,
        body_a: &RigidBody,
        body_b: &RigidBody,
        restitution: f32,
        dt: f32,
    ) -> Self {
        let point = (contact.point_a + contact.point_b) * 0.5;
        let normal = contact.normal;

        // Compute tangent directions
        let (tangent1, tangent2) = compute_tangent_basis(normal);

        // Radius vectors
        let r_a = point - body_a.position;
        let r_b = point - body_b.position;

        // Compute effective masses
        let normal_mass = compute_effective_mass(body_a, body_b, r_a, r_b, normal);
        let tangent1_mass = compute_effective_mass(body_a, body_b, r_a, r_b, tangent1);
        let tangent2_mass = compute_effective_mass(body_a, body_b, r_a, r_b, tangent2);

        // Compute velocity bias for restitution
        let relative_velocity = compute_relative_velocity(body_a, body_b, r_a, r_b);
        let normal_velocity = relative_velocity.dot(normal);

        // Restitution velocity bias - only apply for significant impacts
        // Use higher threshold to avoid jitter on slow contacts
        let restitution_threshold = 2.0; // m/s - only bounce for faster impacts
        let restitution_bias = if normal_velocity < -restitution_threshold {
            -restitution * normal_velocity
        } else {
            0.0
        };

        // Baumgarte stabilization - soft constraint to push objects apart
        // Use moderate bias scaled by dt to correct penetration gradually
        let penetration = contact.depth;
        let slop = 0.01;
        let erp = 0.1; // Error reduction - 10% correction per frame
        let max_bias = 5.0; // Cap to prevent explosion

        let penetration_bias = if penetration > slop {
            let bias = erp * (penetration - slop) / dt;
            bias.min(max_bias)
        } else {
            0.0
        };

        let velocity_bias = restitution_bias + penetration_bias;

        Self {
            point,
            normal,
            tangent1,
            tangent2,
            r_a,
            r_b,
            normal_mass,
            tangent1_mass,
            tangent2_mass,
            velocity_bias,
            normal_impulse: contact.normal_impulse,
            tangent1_impulse: contact.tangent_impulse_1,
            tangent2_impulse: contact.tangent_impulse_2,
        }
    }

    /// Solves the normal constraint (non-penetration)
    pub fn solve_normal(&mut self, body_a: &mut RigidBody, body_b: &mut RigidBody) {
        // Compute relative velocity at contact point
        let relative_velocity = compute_relative_velocity(body_a, body_b, self.r_a, self.r_b);
        let normal_velocity = relative_velocity.dot(self.normal);

        // Compute impulse
        let mut impulse = self.normal_mass * (-normal_velocity + self.velocity_bias);

        // Clamp accumulated impulse (non-negative for contacts)
        let old_impulse = self.normal_impulse;
        self.normal_impulse = (old_impulse + impulse).max(0.0);
        impulse = self.normal_impulse - old_impulse;

        // Apply impulse
        let p = self.normal * impulse;
        apply_impulse(body_a, body_b, p, self.r_a, self.r_b);
    }

    /// Solves the friction constraints
    pub fn solve_friction(&mut self, body_a: &mut RigidBody, body_b: &mut RigidBody, friction: f32) {
        let max_friction = friction * self.normal_impulse;

        // Tangent 1
        {
            let relative_velocity = compute_relative_velocity(body_a, body_b, self.r_a, self.r_b);
            let tangent_velocity = relative_velocity.dot(self.tangent1);

            let mut impulse = self.tangent1_mass * (-tangent_velocity);

            // Clamp accumulated impulse
            let old_impulse = self.tangent1_impulse;
            self.tangent1_impulse = (old_impulse + impulse).clamp(-max_friction, max_friction);
            impulse = self.tangent1_impulse - old_impulse;

            let p = self.tangent1 * impulse;
            apply_impulse(body_a, body_b, p, self.r_a, self.r_b);
        }

        // Tangent 2
        {
            let relative_velocity = compute_relative_velocity(body_a, body_b, self.r_a, self.r_b);
            let tangent_velocity = relative_velocity.dot(self.tangent2);

            let mut impulse = self.tangent2_mass * (-tangent_velocity);

            // Clamp accumulated impulse
            let old_impulse = self.tangent2_impulse;
            self.tangent2_impulse = (old_impulse + impulse).clamp(-max_friction, max_friction);
            impulse = self.tangent2_impulse - old_impulse;

            let p = self.tangent2 * impulse;
            apply_impulse(body_a, body_b, p, self.r_a, self.r_b);
        }
    }

    /// Stores accumulated impulses back to the contact point
    pub fn store_impulses(&self, contact: &mut ContactPoint) {
        contact.normal_impulse = self.normal_impulse;
        contact.tangent_impulse_1 = self.tangent1_impulse;
        contact.tangent_impulse_2 = self.tangent2_impulse;
    }
}

/// Computes an orthonormal tangent basis from a normal
fn compute_tangent_basis(normal: Vec3) -> (Vec3, Vec3) {
    let tangent1 = if normal.x.abs() >= 0.57735 {
        Vec3::new(normal.y, -normal.x, 0.0).normalize()
    } else {
        Vec3::new(0.0, normal.z, -normal.y).normalize()
    };

    let tangent2 = normal.cross(tangent1);

    (tangent1, tangent2)
}

/// Computes the effective mass for a constraint direction
fn compute_effective_mass(
    body_a: &RigidBody,
    body_b: &RigidBody,
    r_a: Vec3,
    r_b: Vec3,
    direction: Vec3,
) -> f32 {
    let rn_a = r_a.cross(direction);
    let rn_b = r_b.cross(direction);

    let k = body_a.inv_mass
        + body_b.inv_mass
        + rn_a.dot(body_a.inv_inertia_world * rn_a)
        + rn_b.dot(body_b.inv_inertia_world * rn_b);

    if k > 0.0 {
        1.0 / k
    } else {
        0.0
    }
}

/// Computes the relative velocity at the contact point
/// Returns vel_b - vel_a (velocity of B relative to A)
/// When normal points from B to A, positive dot means separating, negative means approaching
fn compute_relative_velocity(body_a: &RigidBody, body_b: &RigidBody, r_a: Vec3, r_b: Vec3) -> Vec3 {
    let vel_a = body_a.linear_velocity + body_a.angular_velocity.cross(r_a);
    let vel_b = body_b.linear_velocity + body_b.angular_velocity.cross(r_b);
    vel_b - vel_a
}

/// Applies an impulse to both bodies
/// With normal pointing from B to A, a positive impulse pushes B away from A
fn apply_impulse(body_a: &mut RigidBody, body_b: &mut RigidBody, impulse: Vec3, r_a: Vec3, r_b: Vec3) {
    // Impulse is in direction of normal (from B to A)
    // A receives negative impulse (pushed away from B, opposite to normal)
    // B receives positive impulse (pushed away from A, along normal)
    body_a.linear_velocity -= impulse * body_a.inv_mass;
    body_a.angular_velocity -= body_a.inv_inertia_world * r_a.cross(impulse);

    body_b.linear_velocity += impulse * body_b.inv_mass;
    body_b.angular_velocity += body_b.inv_inertia_world * r_b.cross(impulse);
}

/// Warm starts the constraint by applying previously accumulated impulses
pub fn warm_start(
    constraint: &ContactVelocityConstraint,
    body_a: &mut RigidBody,
    body_b: &mut RigidBody,
) {
    let p = constraint.normal * constraint.normal_impulse
        + constraint.tangent1 * constraint.tangent1_impulse
        + constraint.tangent2 * constraint.tangent2_impulse;

    apply_impulse(body_a, body_b, p, constraint.r_a, constraint.r_b);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dynamics::BodyType;

    #[test]
    fn test_tangent_basis() {
        let normal = Vec3::Y;
        let (t1, t2) = compute_tangent_basis(normal);

        // Tangents should be perpendicular to normal
        assert!(t1.dot(normal).abs() < 0.0001);
        assert!(t2.dot(normal).abs() < 0.0001);

        // Tangents should be perpendicular to each other
        assert!(t1.dot(t2).abs() < 0.0001);

        // Tangents should be unit vectors
        assert!((t1.length() - 1.0).abs() < 0.0001);
        assert!((t2.length() - 1.0).abs() < 0.0001);
    }

    #[test]
    fn test_effective_mass() {
        let mut body_a = RigidBody::new().with_mass(1.0);
        let body_b = RigidBody::new().with_type(BodyType::Static);

        body_a.update_world_inertia();

        let r_a = Vec3::ZERO;
        let r_b = Vec3::ZERO;
        let direction = Vec3::Y;

        let mass = compute_effective_mass(&body_a, &body_b, r_a, r_b, direction);

        // With body_b static and r = 0, effective mass should equal body_a's mass
        assert!((mass - 1.0).abs() < 0.1);
    }
}
