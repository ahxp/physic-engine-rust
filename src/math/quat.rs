use std::ops::{Mul, MulAssign, Neg};

use super::vec3::Vec3;

/// A quaternion representing a rotation in 3D space.
///
/// Stored as (x, y, z, w) where w is the scalar part.
/// Always kept normalized for rotation operations.
#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(C)]
pub struct Quat {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub w: f32,
}

impl Default for Quat {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Quat {
    /// Identity quaternion (no rotation)
    pub const IDENTITY: Self = Self::new(0.0, 0.0, 0.0, 1.0);

    /// Creates a new quaternion from components
    #[inline]
    pub const fn new(x: f32, y: f32, z: f32, w: f32) -> Self {
        Self { x, y, z, w }
    }

    /// Creates a quaternion from a rotation axis and angle (in radians)
    #[inline]
    pub fn from_axis_angle(axis: Vec3, angle: f32) -> Self {
        let half_angle = angle * 0.5;
        let s = half_angle.sin();
        let c = half_angle.cos();
        let axis = axis.normalize();
        Self::new(axis.x * s, axis.y * s, axis.z * s, c)
    }

    /// Creates a quaternion from Euler angles (in radians)
    /// Order: X (pitch) -> Y (yaw) -> Z (roll)
    #[inline]
    pub fn from_euler(pitch: f32, yaw: f32, roll: f32) -> Self {
        let (sp, cp) = (pitch * 0.5).sin_cos();
        let (sy, cy) = (yaw * 0.5).sin_cos();
        let (sr, cr) = (roll * 0.5).sin_cos();

        Self::new(
            sp * cy * cr - cp * sy * sr,
            cp * sy * cr + sp * cy * sr,
            cp * cy * sr - sp * sy * cr,
            cp * cy * cr + sp * sy * sr,
        )
    }

    /// Creates a quaternion that rotates from one vector to another
    #[inline]
    pub fn from_rotation_arc(from: Vec3, to: Vec3) -> Self {
        let from = from.normalize();
        let to = to.normalize();

        let dot = from.dot(to);

        if dot > 0.9999 {
            // Vectors are nearly parallel
            return Self::IDENTITY;
        }

        if dot < -0.9999 {
            // Vectors are nearly opposite
            let axis = Vec3::X.cross(from);
            let axis = if axis.length_squared() < 0.0001 {
                Vec3::Y.cross(from)
            } else {
                axis
            };
            return Self::from_axis_angle(axis.normalize(), std::f32::consts::PI);
        }

        let cross = from.cross(to);
        Self::new(cross.x, cross.y, cross.z, 1.0 + dot).normalize()
    }

    /// Returns the rotation axis and angle (in radians)
    #[inline]
    pub fn to_axis_angle(self) -> (Vec3, f32) {
        let angle = 2.0 * self.w.clamp(-1.0, 1.0).acos();
        let s = (1.0 - self.w * self.w).sqrt();

        if s < 0.0001 {
            (Vec3::X, angle)
        } else {
            let axis = Vec3::new(self.x / s, self.y / s, self.z / s);
            (axis, angle)
        }
    }

    /// Converts to Euler angles (pitch, yaw, roll) in radians
    #[inline]
    pub fn to_euler(self) -> (f32, f32, f32) {
        // Pitch (x-axis rotation)
        let sinp = 2.0 * (self.w * self.x + self.y * self.z);
        let cosp = 1.0 - 2.0 * (self.x * self.x + self.y * self.y);
        let pitch = sinp.atan2(cosp);

        // Yaw (y-axis rotation)
        let siny = 2.0 * (self.w * self.y - self.z * self.x);
        let yaw = if siny.abs() >= 1.0 {
            std::f32::consts::FRAC_PI_2.copysign(siny)
        } else {
            siny.asin()
        };

        // Roll (z-axis rotation)
        let sinr = 2.0 * (self.w * self.z + self.x * self.y);
        let cosr = 1.0 - 2.0 * (self.y * self.y + self.z * self.z);
        let roll = sinr.atan2(cosr);

        (pitch, yaw, roll)
    }

    /// Returns the squared length of the quaternion
    #[inline]
    pub fn length_squared(self) -> f32 {
        self.x * self.x + self.y * self.y + self.z * self.z + self.w * self.w
    }

    /// Returns the length of the quaternion
    #[inline]
    pub fn length(self) -> f32 {
        self.length_squared().sqrt()
    }

    /// Returns a normalized quaternion
    #[inline]
    pub fn normalize(self) -> Self {
        let len = self.length();
        if len > 1e-10 {
            let inv_len = 1.0 / len;
            Self::new(
                self.x * inv_len,
                self.y * inv_len,
                self.z * inv_len,
                self.w * inv_len,
            )
        } else {
            Self::IDENTITY
        }
    }

    /// Returns the conjugate (inverse rotation for unit quaternions)
    #[inline]
    pub fn conjugate(self) -> Self {
        Self::new(-self.x, -self.y, -self.z, self.w)
    }

    /// Returns the inverse of the quaternion
    #[inline]
    pub fn inverse(self) -> Self {
        let len_sq = self.length_squared();
        if len_sq > 1e-10 {
            let inv_len_sq = 1.0 / len_sq;
            Self::new(
                -self.x * inv_len_sq,
                -self.y * inv_len_sq,
                -self.z * inv_len_sq,
                self.w * inv_len_sq,
            )
        } else {
            Self::IDENTITY
        }
    }

    /// Dot product of two quaternions
    #[inline]
    pub fn dot(self, other: Self) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z + self.w * other.w
    }

    /// Rotates a vector by this quaternion
    #[inline]
    pub fn rotate_vec(self, v: Vec3) -> Vec3 {
        // Optimized quaternion-vector rotation
        let qv = Vec3::new(self.x, self.y, self.z);
        let uv = qv.cross(v);
        let uuv = qv.cross(uv);
        v + (uv * self.w + uuv) * 2.0
    }

    /// Inverse rotates a vector (rotates by conjugate)
    #[inline]
    pub fn inverse_rotate_vec(self, v: Vec3) -> Vec3 {
        self.conjugate().rotate_vec(v)
    }

    /// Spherical linear interpolation between two quaternions
    #[inline]
    pub fn slerp(self, other: Self, t: f32) -> Self {
        let mut dot = self.dot(other);

        // If dot < 0, negate one quaternion to take shorter path
        let other = if dot < 0.0 {
            dot = -dot;
            -other
        } else {
            other
        };

        // If quaternions are very close, use linear interpolation
        if dot > 0.9995 {
            return Self::new(
                self.x + t * (other.x - self.x),
                self.y + t * (other.y - self.y),
                self.z + t * (other.z - self.z),
                self.w + t * (other.w - self.w),
            )
            .normalize();
        }

        let theta = dot.clamp(-1.0, 1.0).acos();
        let theta_t = theta * t;

        let sin_theta = theta.sin();
        let sin_theta_t = theta_t.sin();
        let cos_theta_t = theta_t.cos();

        let s0 = cos_theta_t - dot * sin_theta_t / sin_theta;
        let s1 = sin_theta_t / sin_theta;

        Self::new(
            s0 * self.x + s1 * other.x,
            s0 * self.y + s1 * other.y,
            s0 * self.z + s1 * other.z,
            s0 * self.w + s1 * other.w,
        )
    }

    /// Normalized linear interpolation (faster but less accurate than slerp)
    #[inline]
    pub fn nlerp(self, other: Self, t: f32) -> Self {
        let mut dot = self.dot(other);
        let other = if dot < 0.0 {
            dot = -dot;
            -other
        } else {
            other
        };

        Self::new(
            self.x + t * (other.x - self.x),
            self.y + t * (other.y - self.y),
            self.z + t * (other.z - self.z),
            self.w + t * (other.w - self.w),
        )
        .normalize()
    }

    /// Returns the local X axis after rotation
    #[inline]
    pub fn local_x(self) -> Vec3 {
        self.rotate_vec(Vec3::X)
    }

    /// Returns the local Y axis after rotation
    #[inline]
    pub fn local_y(self) -> Vec3 {
        self.rotate_vec(Vec3::Y)
    }

    /// Returns the local Z axis after rotation
    #[inline]
    pub fn local_z(self) -> Vec3 {
        self.rotate_vec(Vec3::Z)
    }

    /// Integrates angular velocity over a time step
    /// Returns the new orientation after rotating by angular_velocity * dt
    #[inline]
    pub fn integrate(self, angular_velocity: Vec3, dt: f32) -> Self {
        let omega = angular_velocity * dt;
        let omega_len = omega.length();

        if omega_len < 1e-10 {
            return self;
        }

        let half_omega_len = omega_len * 0.5;
        let s = half_omega_len.sin() / omega_len;
        let c = half_omega_len.cos();

        let dq = Quat::new(omega.x * s, omega.y * s, omega.z * s, c);

        (dq * self).normalize()
    }
}

impl Mul for Quat {
    type Output = Self;

    /// Quaternion multiplication (combines rotations)
    #[inline]
    fn mul(self, other: Self) -> Self {
        Self::new(
            self.w * other.x + self.x * other.w + self.y * other.z - self.z * other.y,
            self.w * other.y - self.x * other.z + self.y * other.w + self.z * other.x,
            self.w * other.z + self.x * other.y - self.y * other.x + self.z * other.w,
            self.w * other.w - self.x * other.x - self.y * other.y - self.z * other.z,
        )
    }
}

impl MulAssign for Quat {
    #[inline]
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}

impl Neg for Quat {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z, -self.w)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    const EPSILON: f32 = 1e-5;

    fn approx_eq(a: f32, b: f32) -> bool {
        (a - b).abs() < EPSILON
    }

    fn quat_approx_eq(a: Quat, b: Quat) -> bool {
        // Quaternions q and -q represent the same rotation
        let dot = a.dot(b);
        dot.abs() > 1.0 - EPSILON
    }

    fn vec3_approx_eq(a: Vec3, b: Vec3) -> bool {
        approx_eq(a.x, b.x) && approx_eq(a.y, b.y) && approx_eq(a.z, b.z)
    }

    #[test]
    fn test_identity() {
        let q = Quat::IDENTITY;
        let v = Vec3::new(1.0, 2.0, 3.0);
        let rotated = q.rotate_vec(v);
        assert!(vec3_approx_eq(rotated, v));
    }

    #[test]
    fn test_axis_angle() {
        // 90 degree rotation around Z axis
        let q = Quat::from_axis_angle(Vec3::Z, PI / 2.0);
        let v = Vec3::X;
        let rotated = q.rotate_vec(v);
        assert!(vec3_approx_eq(rotated, Vec3::Y));
    }

    #[test]
    fn test_inverse() {
        let q = Quat::from_axis_angle(Vec3::new(1.0, 1.0, 1.0).normalize(), PI / 3.0);
        let v = Vec3::new(1.0, 2.0, 3.0);

        let rotated = q.rotate_vec(v);
        let back = q.inverse().rotate_vec(rotated);

        assert!(vec3_approx_eq(back, v));
    }

    #[test]
    fn test_multiplication() {
        // Two 90 degree rotations around Z should equal one 180 degree rotation
        let q1 = Quat::from_axis_angle(Vec3::Z, PI / 2.0);
        let q2 = q1 * q1;
        let q180 = Quat::from_axis_angle(Vec3::Z, PI);

        assert!(quat_approx_eq(q2, q180));
    }

    #[test]
    fn test_slerp() {
        let q1 = Quat::IDENTITY;
        let q2 = Quat::from_axis_angle(Vec3::Z, PI / 2.0);

        let mid = q1.slerp(q2, 0.5);
        let expected = Quat::from_axis_angle(Vec3::Z, PI / 4.0);

        assert!(quat_approx_eq(mid, expected));
    }

    #[test]
    fn test_normalize() {
        let q = Quat::new(1.0, 2.0, 3.0, 4.0);
        let n = q.normalize();
        assert!(approx_eq(n.length(), 1.0));
    }

    #[test]
    fn test_rotate_vec() {
        // 180 degree rotation around X should flip Y and Z signs
        let q = Quat::from_axis_angle(Vec3::X, PI);
        let v = Vec3::new(0.0, 1.0, 1.0);
        let rotated = q.rotate_vec(v);
        assert!(vec3_approx_eq(rotated, Vec3::new(0.0, -1.0, -1.0)));
    }

    #[test]
    fn test_local_axes() {
        let q = Quat::from_axis_angle(Vec3::Z, PI / 2.0);

        // After 90 deg rotation around Z:
        // Local X should point along world Y
        // Local Y should point along world -X
        assert!(vec3_approx_eq(q.local_x(), Vec3::Y));
        assert!(vec3_approx_eq(q.local_y(), -Vec3::X));
        assert!(vec3_approx_eq(q.local_z(), Vec3::Z));
    }

    #[test]
    fn test_integrate() {
        let q = Quat::IDENTITY;
        let omega = Vec3::new(0.0, 0.0, PI); // PI rad/s around Z
        let dt = 0.5; // After 0.5s, should rotate PI/2 radians

        let result = q.integrate(omega, dt);
        let expected = Quat::from_axis_angle(Vec3::Z, PI / 2.0);

        assert!(quat_approx_eq(result, expected));
    }
}
