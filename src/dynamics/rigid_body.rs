use crate::collision::BodyHandle;
use crate::geometry::{MassProperties, Shape};
use crate::math::{Mat3, Quat, Transform, Vec3};

/// The type of rigid body
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BodyType {
    /// Dynamic bodies are affected by forces and collisions
    Dynamic,
    /// Static bodies never move
    Static,
    /// Kinematic bodies move according to their velocity but aren't affected by forces
    Kinematic,
}

impl Default for BodyType {
    fn default() -> Self {
        BodyType::Dynamic
    }
}

/// A rigid body in the physics simulation
#[derive(Debug, Clone)]
pub struct RigidBody {
    /// Body handle for identification
    pub handle: BodyHandle,
    /// Body type (dynamic, static, kinematic)
    pub body_type: BodyType,

    // Transform
    /// Position in world space
    pub position: Vec3,
    /// Rotation as quaternion
    pub rotation: Quat,

    // Velocities
    /// Linear velocity
    pub linear_velocity: Vec3,
    /// Angular velocity (in radians per second)
    pub angular_velocity: Vec3,

    // Mass properties
    /// Inverse mass (0 for infinite mass / static)
    pub inv_mass: f32,
    /// Local inverse inertia tensor
    pub inv_inertia_local: Mat3,
    /// World space inverse inertia tensor (updated each frame)
    pub inv_inertia_world: Mat3,

    // Forces
    /// Accumulated force (reset each step)
    pub force: Vec3,
    /// Accumulated torque (reset each step)
    pub torque: Vec3,

    // Material properties
    /// Friction coefficient
    pub friction: f32,
    /// Restitution (bounciness)
    pub restitution: f32,

    // Damping
    /// Linear damping (0-1)
    pub linear_damping: f32,
    /// Angular damping (0-1)
    pub angular_damping: f32,

    // Collision shape (optional)
    /// The collision shape attached to this body
    pub shape: Option<Shape>,

    // Flags
    /// Whether the body is awake
    pub is_awake: bool,
    /// Whether the body can sleep
    pub can_sleep: bool,
    /// Sleep counter
    sleep_time: f32,

    // User data
    /// Optional user data
    pub user_data: u64,
}

impl Default for RigidBody {
    fn default() -> Self {
        Self {
            handle: BodyHandle::INVALID,
            body_type: BodyType::Dynamic,
            position: Vec3::ZERO,
            rotation: Quat::IDENTITY,
            linear_velocity: Vec3::ZERO,
            angular_velocity: Vec3::ZERO,
            inv_mass: 1.0,
            inv_inertia_local: Mat3::IDENTITY,
            inv_inertia_world: Mat3::IDENTITY,
            force: Vec3::ZERO,
            torque: Vec3::ZERO,
            friction: 0.6,         // Good friction for stability
            restitution: 0.3,      // Moderate bounce
            linear_damping: 0.0,   // No artificial linear damping
            angular_damping: 0.1,  // Strong angular damping for stability
            shape: None,
            is_awake: true,
            can_sleep: true,
            sleep_time: 0.0,
            user_data: 0,
        }
    }
}

impl RigidBody {
    /// Creates a new rigid body
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the body type
    pub fn with_type(mut self, body_type: BodyType) -> Self {
        self.body_type = body_type;
        if body_type != BodyType::Dynamic {
            self.inv_mass = 0.0;
            self.inv_inertia_local = Mat3::ZERO;
            self.inv_inertia_world = Mat3::ZERO;
        }
        self
    }

    /// Sets the position
    pub fn with_position(mut self, position: Vec3) -> Self {
        self.position = position;
        self
    }

    /// Sets the rotation
    pub fn with_rotation(mut self, rotation: Quat) -> Self {
        self.rotation = rotation;
        self.update_world_inertia();
        self
    }

    /// Sets the mass (automatically computes inverse mass)
    pub fn with_mass(mut self, mass: f32) -> Self {
        if mass > 0.0 && self.body_type == BodyType::Dynamic {
            self.inv_mass = 1.0 / mass;
        } else {
            self.inv_mass = 0.0;
        }
        self
    }

    /// Sets mass properties from a shape and density
    pub fn with_mass_properties(mut self, props: MassProperties) -> Self {
        if self.body_type != BodyType::Dynamic {
            return self;
        }

        self.inv_mass = props.inv_mass();
        self.inv_inertia_local = props.inv_inertia();
        self.update_world_inertia();
        self
    }

    /// Sets the collision shape
    pub fn with_shape(mut self, shape: Shape) -> Self {
        self.shape = Some(shape);
        self
    }

    /// Sets friction
    pub fn with_friction(mut self, friction: f32) -> Self {
        self.friction = friction.clamp(0.0, 1.0);
        self
    }

    /// Sets restitution
    pub fn with_restitution(mut self, restitution: f32) -> Self {
        self.restitution = restitution.clamp(0.0, 1.0);
        self
    }

    /// Sets linear damping
    pub fn with_linear_damping(mut self, damping: f32) -> Self {
        self.linear_damping = damping.clamp(0.0, 1.0);
        self
    }

    /// Sets angular damping
    pub fn with_angular_damping(mut self, damping: f32) -> Self {
        self.angular_damping = damping.clamp(0.0, 1.0);
        self
    }

    /// Returns the transform of this body
    pub fn transform(&self) -> Transform {
        Transform::new(self.position, self.rotation)
    }

    /// Returns the mass (inverse of inv_mass, or infinity for static)
    pub fn mass(&self) -> f32 {
        if self.inv_mass > 0.0 {
            1.0 / self.inv_mass
        } else {
            f32::INFINITY
        }
    }

    /// Returns true if this body has finite mass
    pub fn has_finite_mass(&self) -> bool {
        self.inv_mass > 0.0
    }

    /// Returns true if this is a dynamic body
    pub fn is_dynamic(&self) -> bool {
        self.body_type == BodyType::Dynamic
    }

    /// Returns true if this is a static body
    pub fn is_static(&self) -> bool {
        self.body_type == BodyType::Static
    }

    /// Returns true if this is a kinematic body
    pub fn is_kinematic(&self) -> bool {
        self.body_type == BodyType::Kinematic
    }

    /// Applies a force at the center of mass
    pub fn apply_force(&mut self, force: Vec3) {
        if self.is_dynamic() {
            self.force += force;
            self.wake_up();
        }
    }

    /// Applies a force at a world point
    pub fn apply_force_at_point(&mut self, force: Vec3, point: Vec3) {
        if self.is_dynamic() {
            self.force += force;
            self.torque += (point - self.position).cross(force);
            self.wake_up();
        }
    }

    /// Applies a torque
    pub fn apply_torque(&mut self, torque: Vec3) {
        if self.is_dynamic() {
            self.torque += torque;
            self.wake_up();
        }
    }

    /// Applies an impulse at the center of mass
    pub fn apply_impulse(&mut self, impulse: Vec3) {
        if self.is_dynamic() {
            self.linear_velocity += impulse * self.inv_mass;
            self.wake_up();
        }
    }

    /// Applies an impulse at a world point
    pub fn apply_impulse_at_point(&mut self, impulse: Vec3, point: Vec3) {
        if self.is_dynamic() {
            self.linear_velocity += impulse * self.inv_mass;
            let r = point - self.position;
            self.angular_velocity += self.inv_inertia_world * r.cross(impulse);
            self.wake_up();
        }
    }

    /// Applies an angular impulse
    pub fn apply_angular_impulse(&mut self, impulse: Vec3) {
        if self.is_dynamic() {
            self.angular_velocity += self.inv_inertia_world * impulse;
            self.wake_up();
        }
    }

    /// Gets the velocity at a world point
    pub fn velocity_at_point(&self, point: Vec3) -> Vec3 {
        self.linear_velocity + self.angular_velocity.cross(point - self.position)
    }

    /// Clears accumulated forces
    pub fn clear_forces(&mut self) {
        self.force = Vec3::ZERO;
        self.torque = Vec3::ZERO;
    }

    /// Updates the world space inertia tensor
    pub fn update_world_inertia(&mut self) {
        if self.body_type != BodyType::Dynamic {
            self.inv_inertia_world = Mat3::ZERO;
            return;
        }

        let rot = Mat3::from_quat(self.rotation);
        self.inv_inertia_world = rot * self.inv_inertia_local * rot.transpose();
    }

    /// Wakes up the body
    pub fn wake_up(&mut self) {
        self.is_awake = true;
        self.sleep_time = 0.0;
    }

    /// Puts the body to sleep
    pub fn sleep(&mut self) {
        if self.can_sleep {
            self.is_awake = false;
            self.linear_velocity = Vec3::ZERO;
            self.angular_velocity = Vec3::ZERO;
            self.force = Vec3::ZERO;
            self.torque = Vec3::ZERO;
        }
    }

    /// Updates sleep state
    pub fn update_sleep(&mut self, dt: f32, sleep_threshold: f32, sleep_time_threshold: f32) {
        if !self.can_sleep || !self.is_dynamic() {
            return;
        }

        let kinetic_energy = self.linear_velocity.length_squared()
            + self.angular_velocity.length_squared();

        if kinetic_energy < sleep_threshold {
            self.sleep_time += dt;
            if self.sleep_time > sleep_time_threshold {
                self.sleep();
            }
        } else {
            self.sleep_time = 0.0;
        }
    }
}

/// Description for creating a rigid body
#[derive(Debug, Clone)]
pub struct RigidBodyDesc {
    pub body_type: BodyType,
    pub position: Vec3,
    pub rotation: Quat,
    pub linear_velocity: Vec3,
    pub angular_velocity: Vec3,
    pub mass: f32,
    pub friction: f32,
    pub restitution: f32,
    pub linear_damping: f32,
    pub angular_damping: f32,
    pub can_sleep: bool,
    pub user_data: u64,
}

impl Default for RigidBodyDesc {
    fn default() -> Self {
        Self {
            body_type: BodyType::Dynamic,
            position: Vec3::ZERO,
            rotation: Quat::IDENTITY,
            linear_velocity: Vec3::ZERO,
            angular_velocity: Vec3::ZERO,
            mass: 1.0,
            friction: 0.6,          // Good friction for stability
            restitution: 0.3,       // Moderate bounce
            linear_damping: 0.0,    // No artificial linear damping
            angular_damping: 0.1,   // Strong angular damping for stability
            can_sleep: true,
            user_data: 0,
        }
    }
}

impl RigidBodyDesc {
    /// Creates a new dynamic body description
    pub fn dynamic() -> Self {
        Self::default()
    }

    /// Creates a new static body description
    pub fn fixed() -> Self {
        Self {
            body_type: BodyType::Static,
            mass: 0.0,
            ..Self::default()
        }
    }

    /// Creates a new kinematic body description
    pub fn kinematic() -> Self {
        Self {
            body_type: BodyType::Kinematic,
            mass: 0.0,
            ..Self::default()
        }
    }

    /// Sets the position
    pub fn with_position(mut self, position: Vec3) -> Self {
        self.position = position;
        self
    }

    /// Sets the rotation
    pub fn with_rotation(mut self, rotation: Quat) -> Self {
        self.rotation = rotation;
        self
    }

    /// Sets the mass
    pub fn with_mass(mut self, mass: f32) -> Self {
        self.mass = mass;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_body_creation() {
        let body = RigidBody::new()
            .with_position(Vec3::new(1.0, 2.0, 3.0))
            .with_mass(2.0);

        assert_eq!(body.position, Vec3::new(1.0, 2.0, 3.0));
        assert!((body.inv_mass - 0.5).abs() < 0.0001);
    }

    #[test]
    fn test_static_body() {
        let body = RigidBody::new().with_type(BodyType::Static);

        assert!(body.is_static());
        assert!(!body.has_finite_mass());
        assert_eq!(body.inv_mass, 0.0);
    }

    #[test]
    fn test_apply_impulse() {
        let mut body = RigidBody::new().with_mass(1.0);

        body.apply_impulse(Vec3::new(1.0, 0.0, 0.0));

        assert_eq!(body.linear_velocity, Vec3::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn test_velocity_at_point() {
        let mut body = RigidBody::new()
            .with_position(Vec3::ZERO)
            .with_mass(1.0);

        body.linear_velocity = Vec3::new(1.0, 0.0, 0.0);
        body.angular_velocity = Vec3::new(0.0, 0.0, 1.0); // Rotating around Z

        // Point at (0, 1, 0) should have velocity (1, 0, 0) + (0, 0, 1) × (0, 1, 0) = (1, 0, 0) + (-1, 0, 0) = (0, 0, 0)
        // Wait, cross product: (0, 0, 1) × (0, 1, 0) = (-1, 0, 0)
        let vel = body.velocity_at_point(Vec3::new(0.0, 1.0, 0.0));
        assert!((vel - Vec3::new(0.0, 0.0, 0.0)).length() < 0.0001);
    }
}
