use super::mat3::Mat3;
use super::quat::Quat;
use super::vec3::Vec3;

/// A rigid body transformation combining position and rotation.
///
/// Represents a coordinate frame in 3D space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Transform {
    /// Position (translation)
    pub position: Vec3,
    /// Rotation (as quaternion)
    pub rotation: Quat,
}

impl Default for Transform {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Transform {
    /// Identity transform (no translation or rotation)
    pub const IDENTITY: Self = Self {
        position: Vec3::ZERO,
        rotation: Quat::IDENTITY,
    };

    /// Creates a new transform from position and rotation
    #[inline]
    pub const fn new(position: Vec3, rotation: Quat) -> Self {
        Self { position, rotation }
    }

    /// Creates a transform with only translation
    #[inline]
    pub const fn from_position(position: Vec3) -> Self {
        Self {
            position,
            rotation: Quat::IDENTITY,
        }
    }

    /// Creates a transform with only rotation
    #[inline]
    pub const fn from_rotation(rotation: Quat) -> Self {
        Self {
            position: Vec3::ZERO,
            rotation,
        }
    }

    /// Creates a transform from a rotation matrix
    #[inline]
    pub fn from_rotation_matrix(position: Vec3, rotation: Mat3) -> Self {
        Self {
            position,
            rotation: rotation.to_quat(),
        }
    }

    /// Returns the rotation as a 3x3 matrix
    #[inline]
    pub fn rotation_matrix(self) -> Mat3 {
        Mat3::from_quat(self.rotation)
    }

    /// Transforms a point from local space to world space
    #[inline]
    pub fn transform_point(self, point: Vec3) -> Vec3 {
        self.rotation.rotate_vec(point) + self.position
    }

    /// Transforms a vector (direction) from local space to world space
    /// Unlike points, vectors are not affected by translation
    #[inline]
    pub fn transform_vector(self, vector: Vec3) -> Vec3 {
        self.rotation.rotate_vec(vector)
    }

    /// Inverse transforms a point from world space to local space
    #[inline]
    pub fn inverse_transform_point(self, point: Vec3) -> Vec3 {
        self.rotation.inverse_rotate_vec(point - self.position)
    }

    /// Inverse transforms a vector from world space to local space
    #[inline]
    pub fn inverse_transform_vector(self, vector: Vec3) -> Vec3 {
        self.rotation.inverse_rotate_vec(vector)
    }

    /// Returns the inverse of this transform
    #[inline]
    pub fn inverse(self) -> Self {
        let inv_rotation = self.rotation.conjugate();
        let inv_position = inv_rotation.rotate_vec(-self.position);
        Self {
            position: inv_position,
            rotation: inv_rotation,
        }
    }

    /// Combines two transforms: self * other
    /// The result transforms from other's local space to self's local space to world
    #[inline]
    pub fn compose(self, other: Self) -> Self {
        Self {
            position: self.transform_point(other.position),
            rotation: self.rotation * other.rotation,
        }
    }

    /// Returns the local X axis in world space
    #[inline]
    pub fn local_x(self) -> Vec3 {
        self.rotation.local_x()
    }

    /// Returns the local Y axis in world space
    #[inline]
    pub fn local_y(self) -> Vec3 {
        self.rotation.local_y()
    }

    /// Returns the local Z axis in world space
    #[inline]
    pub fn local_z(self) -> Vec3 {
        self.rotation.local_z()
    }

    /// Linear interpolation between two transforms
    #[inline]
    pub fn lerp(self, other: Self, t: f32) -> Self {
        Self {
            position: self.position.lerp(other.position, t),
            rotation: self.rotation.slerp(other.rotation, t),
        }
    }

    /// Creates a transform that looks from `eye` towards `target` with the given `up` vector
    #[inline]
    pub fn look_at(eye: Vec3, target: Vec3, up: Vec3) -> Self {
        let forward = (target - eye).normalize();
        let right = up.cross(forward).normalize();
        let up = forward.cross(right);

        let rotation_matrix = Mat3::from_cols(right, up, forward);
        Self {
            position: eye,
            rotation: rotation_matrix.to_quat(),
        }
    }

    /// Returns true if this transform is approximately equal to another
    #[inline]
    pub fn approx_eq(self, other: Self, epsilon: f32) -> bool {
        self.position.distance_squared(other.position) < epsilon * epsilon
            && self.rotation.dot(other.rotation).abs() > 1.0 - epsilon
    }
}

/// Isometry is an alias for Transform (both represent rigid transformations)
pub type Isometry = Transform;

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    const EPSILON: f32 = 1e-5;

    fn vec3_approx_eq(a: Vec3, b: Vec3) -> bool {
        a.distance_squared(b) < EPSILON * EPSILON
    }

    #[test]
    fn test_identity() {
        let t = Transform::IDENTITY;
        let p = Vec3::new(1.0, 2.0, 3.0);

        assert!(vec3_approx_eq(t.transform_point(p), p));
        assert!(vec3_approx_eq(t.transform_vector(p), p));
    }

    #[test]
    fn test_translation_only() {
        let t = Transform::from_position(Vec3::new(1.0, 2.0, 3.0));
        let p = Vec3::new(1.0, 1.0, 1.0);

        // Points are affected by translation
        assert!(vec3_approx_eq(
            t.transform_point(p),
            Vec3::new(2.0, 3.0, 4.0)
        ));

        // Vectors are not affected by translation
        assert!(vec3_approx_eq(t.transform_vector(p), p));
    }

    #[test]
    fn test_rotation_only() {
        // 90 degree rotation around Z
        let t = Transform::from_rotation(Quat::from_axis_angle(Vec3::Z, PI / 2.0));

        let p = Vec3::X;
        let rotated = t.transform_point(p);
        assert!(vec3_approx_eq(rotated, Vec3::Y));
    }

    #[test]
    fn test_combined() {
        let t = Transform::new(
            Vec3::new(1.0, 0.0, 0.0),
            Quat::from_axis_angle(Vec3::Z, PI / 2.0),
        );

        // First rotate X -> Y, then translate by (1, 0, 0)
        let p = Vec3::X;
        let result = t.transform_point(p);
        assert!(vec3_approx_eq(result, Vec3::new(1.0, 1.0, 0.0)));
    }

    #[test]
    fn test_inverse() {
        let t = Transform::new(
            Vec3::new(1.0, 2.0, 3.0),
            Quat::from_axis_angle(Vec3::new(1.0, 1.0, 1.0).normalize(), PI / 4.0),
        );

        let p = Vec3::new(4.0, 5.0, 6.0);
        let transformed = t.transform_point(p);
        let back = t.inverse().transform_point(transformed);

        assert!(vec3_approx_eq(back, p));
    }

    #[test]
    fn test_compose() {
        let t1 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let t2 = Transform::from_rotation(Quat::from_axis_angle(Vec3::Z, PI / 2.0));

        // Compose: first apply t2 (rotation), then t1 (translation)
        let composed = t1.compose(t2);

        let p = Vec3::X;
        // After t2: X -> Y
        // After t1: Y + (1, 0, 0) = (1, 1, 0)
        let result = composed.transform_point(p);
        assert!(vec3_approx_eq(result, Vec3::new(1.0, 1.0, 0.0)));
    }

    #[test]
    fn test_local_axes() {
        let t = Transform::from_rotation(Quat::from_axis_angle(Vec3::Z, PI / 2.0));

        assert!(vec3_approx_eq(t.local_x(), Vec3::Y));
        assert!(vec3_approx_eq(t.local_y(), -Vec3::X));
        assert!(vec3_approx_eq(t.local_z(), Vec3::Z));
    }

    #[test]
    fn test_lerp() {
        let t1 = Transform::IDENTITY;
        let t2 = Transform::new(
            Vec3::new(10.0, 0.0, 0.0),
            Quat::from_axis_angle(Vec3::Z, PI / 2.0), // 90 degrees
        );

        let mid = t1.lerp(t2, 0.5);

        assert!(vec3_approx_eq(mid.position, Vec3::new(5.0, 0.0, 0.0)));

        // Mid rotation should be approximately 45 degrees around Z
        // After 45 degrees rotation around Z, local X points towards (+1, +1, 0) normalized
        let local_x = mid.rotation.local_x();
        let expected = Vec3::new(1.0, 1.0, 0.0).normalize();
        assert!(
            vec3_approx_eq(local_x, expected),
            "local_x = {:?}, expected = {:?}",
            local_x,
            expected
        );
    }
}
