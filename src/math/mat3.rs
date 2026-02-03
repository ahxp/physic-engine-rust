use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

use super::quat::Quat;
use super::vec3::Vec3;

/// A 3x3 matrix stored in column-major order.
///
/// Used for rotation matrices and inertia tensors.
#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(C)]
pub struct Mat3 {
    /// Columns of the matrix
    pub cols: [Vec3; 3],
}

impl Default for Mat3 {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Mat3 {
    /// Zero matrix
    pub const ZERO: Self = Self {
        cols: [Vec3::ZERO, Vec3::ZERO, Vec3::ZERO],
    };

    /// Identity matrix
    pub const IDENTITY: Self = Self {
        cols: [Vec3::X, Vec3::Y, Vec3::Z],
    };

    /// Creates a matrix from column vectors
    #[inline]
    pub const fn from_cols(c0: Vec3, c1: Vec3, c2: Vec3) -> Self {
        Self { cols: [c0, c1, c2] }
    }

    /// Creates a matrix from row vectors
    #[inline]
    pub fn from_rows(r0: Vec3, r1: Vec3, r2: Vec3) -> Self {
        Self::from_cols(
            Vec3::new(r0.x, r1.x, r2.x),
            Vec3::new(r0.y, r1.y, r2.y),
            Vec3::new(r0.z, r1.z, r2.z),
        )
    }

    /// Creates a diagonal matrix
    #[inline]
    pub fn from_diagonal(diag: Vec3) -> Self {
        Self::from_cols(
            Vec3::new(diag.x, 0.0, 0.0),
            Vec3::new(0.0, diag.y, 0.0),
            Vec3::new(0.0, 0.0, diag.z),
        )
    }

    /// Creates a matrix from a scale vector
    #[inline]
    pub fn from_scale(scale: Vec3) -> Self {
        Self::from_diagonal(scale)
    }

    /// Creates a rotation matrix from a quaternion
    #[inline]
    pub fn from_quat(q: Quat) -> Self {
        let x2 = q.x + q.x;
        let y2 = q.y + q.y;
        let z2 = q.z + q.z;

        let xx = q.x * x2;
        let xy = q.x * y2;
        let xz = q.x * z2;
        let yy = q.y * y2;
        let yz = q.y * z2;
        let zz = q.z * z2;
        let wx = q.w * x2;
        let wy = q.w * y2;
        let wz = q.w * z2;

        Self::from_cols(
            Vec3::new(1.0 - (yy + zz), xy + wz, xz - wy),
            Vec3::new(xy - wz, 1.0 - (xx + zz), yz + wx),
            Vec3::new(xz + wy, yz - wx, 1.0 - (xx + yy)),
        )
    }

    /// Creates a rotation matrix from an axis and angle (radians)
    #[inline]
    pub fn from_axis_angle(axis: Vec3, angle: f32) -> Self {
        Self::from_quat(Quat::from_axis_angle(axis, angle))
    }

    /// Converts to a quaternion (assumes this is a rotation matrix)
    #[inline]
    pub fn to_quat(self) -> Quat {
        let trace = self.cols[0].x + self.cols[1].y + self.cols[2].z;

        if trace > 0.0 {
            let s = (trace + 1.0).sqrt() * 2.0;
            let inv_s = 1.0 / s;
            Quat::new(
                (self.cols[1].z - self.cols[2].y) * inv_s,
                (self.cols[2].x - self.cols[0].z) * inv_s,
                (self.cols[0].y - self.cols[1].x) * inv_s,
                s * 0.25,
            )
        } else if self.cols[0].x > self.cols[1].y && self.cols[0].x > self.cols[2].z {
            let s = (1.0 + self.cols[0].x - self.cols[1].y - self.cols[2].z).sqrt() * 2.0;
            let inv_s = 1.0 / s;
            Quat::new(
                s * 0.25,
                (self.cols[0].y + self.cols[1].x) * inv_s,
                (self.cols[2].x + self.cols[0].z) * inv_s,
                (self.cols[1].z - self.cols[2].y) * inv_s,
            )
        } else if self.cols[1].y > self.cols[2].z {
            let s = (1.0 + self.cols[1].y - self.cols[0].x - self.cols[2].z).sqrt() * 2.0;
            let inv_s = 1.0 / s;
            Quat::new(
                (self.cols[0].y + self.cols[1].x) * inv_s,
                s * 0.25,
                (self.cols[1].z + self.cols[2].y) * inv_s,
                (self.cols[2].x - self.cols[0].z) * inv_s,
            )
        } else {
            let s = (1.0 + self.cols[2].z - self.cols[0].x - self.cols[1].y).sqrt() * 2.0;
            let inv_s = 1.0 / s;
            Quat::new(
                (self.cols[2].x + self.cols[0].z) * inv_s,
                (self.cols[1].z + self.cols[2].y) * inv_s,
                s * 0.25,
                (self.cols[0].y - self.cols[1].x) * inv_s,
            )
        }
    }

    /// Returns the transpose of the matrix
    #[inline]
    pub fn transpose(self) -> Self {
        Self::from_cols(
            Vec3::new(self.cols[0].x, self.cols[1].x, self.cols[2].x),
            Vec3::new(self.cols[0].y, self.cols[1].y, self.cols[2].y),
            Vec3::new(self.cols[0].z, self.cols[1].z, self.cols[2].z),
        )
    }

    /// Returns the determinant of the matrix
    #[inline]
    pub fn determinant(self) -> f32 {
        self.cols[0].dot(self.cols[1].cross(self.cols[2]))
    }

    /// Returns the inverse of the matrix, or None if not invertible
    #[inline]
    pub fn try_inverse(self) -> Option<Self> {
        let det = self.determinant();
        if det.abs() < 1e-10 {
            return None;
        }

        let inv_det = 1.0 / det;

        let c0 = self.cols[1].cross(self.cols[2]) * inv_det;
        let c1 = self.cols[2].cross(self.cols[0]) * inv_det;
        let c2 = self.cols[0].cross(self.cols[1]) * inv_det;

        Some(Self::from_cols(
            Vec3::new(c0.x, c1.x, c2.x),
            Vec3::new(c0.y, c1.y, c2.y),
            Vec3::new(c0.z, c1.z, c2.z),
        ))
    }

    /// Returns the inverse of the matrix (panics if not invertible)
    #[inline]
    pub fn inverse(self) -> Self {
        self.try_inverse().expect("Matrix is not invertible")
    }

    /// Transforms a vector by this matrix
    #[inline]
    pub fn transform_vec(self, v: Vec3) -> Vec3 {
        self.cols[0] * v.x + self.cols[1] * v.y + self.cols[2] * v.z
    }

    /// Returns a column of the matrix
    #[inline]
    pub fn col(self, index: usize) -> Vec3 {
        self.cols[index]
    }

    /// Returns a row of the matrix
    #[inline]
    pub fn row(self, index: usize) -> Vec3 {
        Vec3::new(self.cols[0][index], self.cols[1][index], self.cols[2][index])
    }

    /// Returns the diagonal elements
    #[inline]
    pub fn diagonal(self) -> Vec3 {
        Vec3::new(self.cols[0].x, self.cols[1].y, self.cols[2].z)
    }

    /// Returns the trace (sum of diagonal elements)
    #[inline]
    pub fn trace(self) -> f32 {
        self.cols[0].x + self.cols[1].y + self.cols[2].z
    }

    /// Creates a skew-symmetric (cross product) matrix from a vector
    /// Such that skew(v) * u = v.cross(u)
    #[inline]
    pub fn skew(v: Vec3) -> Self {
        Self::from_cols(
            Vec3::new(0.0, v.z, -v.y),
            Vec3::new(-v.z, 0.0, v.x),
            Vec3::new(v.y, -v.x, 0.0),
        )
    }

    /// Creates an outer product matrix (v * w^T)
    #[inline]
    pub fn outer_product(v: Vec3, w: Vec3) -> Self {
        Self::from_cols(v * w.x, v * w.y, v * w.z)
    }

    /// Component-wise multiplication
    #[inline]
    pub fn component_mul(self, other: Self) -> Self {
        Self::from_cols(
            self.cols[0].component_mul(other.cols[0]),
            self.cols[1].component_mul(other.cols[1]),
            self.cols[2].component_mul(other.cols[2]),
        )
    }

    /// Scalar multiplication
    #[inline]
    pub fn scale(self, s: f32) -> Self {
        Self::from_cols(self.cols[0] * s, self.cols[1] * s, self.cols[2] * s)
    }

    /// Returns true if this is approximately equal to another matrix
    #[inline]
    pub fn approx_eq(self, other: Self, epsilon: f32) -> bool {
        (self.cols[0] - other.cols[0]).length_squared() < epsilon * epsilon
            && (self.cols[1] - other.cols[1]).length_squared() < epsilon * epsilon
            && (self.cols[2] - other.cols[2]).length_squared() < epsilon * epsilon
    }
}

impl Add for Mat3 {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self {
        Self::from_cols(
            self.cols[0] + other.cols[0],
            self.cols[1] + other.cols[1],
            self.cols[2] + other.cols[2],
        )
    }
}

impl AddAssign for Mat3 {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        self.cols[0] += other.cols[0];
        self.cols[1] += other.cols[1];
        self.cols[2] += other.cols[2];
    }
}

impl Sub for Mat3 {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self {
        Self::from_cols(
            self.cols[0] - other.cols[0],
            self.cols[1] - other.cols[1],
            self.cols[2] - other.cols[2],
        )
    }
}

impl SubAssign for Mat3 {
    #[inline]
    fn sub_assign(&mut self, other: Self) {
        self.cols[0] -= other.cols[0];
        self.cols[1] -= other.cols[1];
        self.cols[2] -= other.cols[2];
    }
}

impl Mul for Mat3 {
    type Output = Self;

    #[inline]
    fn mul(self, other: Self) -> Self {
        Self::from_cols(
            self.transform_vec(other.cols[0]),
            self.transform_vec(other.cols[1]),
            self.transform_vec(other.cols[2]),
        )
    }
}

impl MulAssign for Mat3 {
    #[inline]
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}

impl Mul<Vec3> for Mat3 {
    type Output = Vec3;

    #[inline]
    fn mul(self, v: Vec3) -> Vec3 {
        self.transform_vec(v)
    }
}

impl Mul<f32> for Mat3 {
    type Output = Self;

    #[inline]
    fn mul(self, s: f32) -> Self {
        self.scale(s)
    }
}

impl Mul<Mat3> for f32 {
    type Output = Mat3;

    #[inline]
    fn mul(self, m: Mat3) -> Mat3 {
        m.scale(self)
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

    fn vec3_approx_eq(a: Vec3, b: Vec3) -> bool {
        approx_eq(a.x, b.x) && approx_eq(a.y, b.y) && approx_eq(a.z, b.z)
    }

    #[test]
    fn test_identity() {
        let m = Mat3::IDENTITY;
        let v = Vec3::new(1.0, 2.0, 3.0);
        assert!(vec3_approx_eq(m * v, v));
    }

    #[test]
    fn test_transpose() {
        let m = Mat3::from_cols(
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(4.0, 5.0, 6.0),
            Vec3::new(7.0, 8.0, 9.0),
        );
        let mt = m.transpose();
        assert!(vec3_approx_eq(mt.col(0), Vec3::new(1.0, 4.0, 7.0)));
        assert!(vec3_approx_eq(mt.col(1), Vec3::new(2.0, 5.0, 8.0)));
        assert!(vec3_approx_eq(mt.col(2), Vec3::new(3.0, 6.0, 9.0)));
    }

    #[test]
    fn test_determinant() {
        let m = Mat3::IDENTITY;
        assert!(approx_eq(m.determinant(), 1.0));

        let m2 = Mat3::from_scale(Vec3::new(2.0, 3.0, 4.0));
        assert!(approx_eq(m2.determinant(), 24.0));
    }

    #[test]
    fn test_inverse() {
        let m = Mat3::from_cols(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 2.0, 0.0),
            Vec3::new(0.0, 0.0, 4.0),
        );
        let inv = m.inverse();
        let identity = m * inv;
        assert!(identity.approx_eq(Mat3::IDENTITY, EPSILON));
    }

    #[test]
    fn test_from_quat() {
        // 90 degree rotation around Z
        let q = Quat::from_axis_angle(Vec3::Z, PI / 2.0);
        let m = Mat3::from_quat(q);

        let v = Vec3::X;
        let rotated_q = q.rotate_vec(v);
        let rotated_m = m * v;

        assert!(vec3_approx_eq(rotated_q, rotated_m));
        assert!(vec3_approx_eq(rotated_m, Vec3::Y));
    }

    #[test]
    fn test_to_quat() {
        let original = Quat::from_axis_angle(Vec3::new(1.0, 1.0, 1.0).normalize(), PI / 3.0);
        let m = Mat3::from_quat(original);
        let recovered = m.to_quat();

        // Compare by rotating a vector
        let v = Vec3::new(1.0, 2.0, 3.0);
        let r1 = original.rotate_vec(v);
        let r2 = recovered.rotate_vec(v);
        assert!(vec3_approx_eq(r1, r2));
    }

    #[test]
    fn test_skew() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let u = Vec3::new(4.0, 5.0, 6.0);

        let skew_v = Mat3::skew(v);
        let result = skew_v * u;
        let expected = v.cross(u);

        assert!(vec3_approx_eq(result, expected));
    }

    #[test]
    fn test_multiplication() {
        let a = Mat3::from_axis_angle(Vec3::Z, PI / 4.0);
        let b = Mat3::from_axis_angle(Vec3::Z, PI / 4.0);
        let c = a * b;

        let expected = Mat3::from_axis_angle(Vec3::Z, PI / 2.0);
        assert!(c.approx_eq(expected, EPSILON));
    }

    #[test]
    fn test_trace() {
        let m = Mat3::from_diagonal(Vec3::new(1.0, 2.0, 3.0));
        assert!(approx_eq(m.trace(), 6.0));
    }
}
