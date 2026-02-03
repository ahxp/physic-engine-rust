use std::ops::{Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign};

/// A 3D vector with f32 components.
///
/// Used throughout the physics engine for positions, velocities, forces, etc.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
#[repr(C)]
pub struct Vec3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Vec3 {
    /// Zero vector (0, 0, 0)
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0);

    /// Unit vector along X axis (1, 0, 0)
    pub const X: Self = Self::new(1.0, 0.0, 0.0);

    /// Unit vector along Y axis (0, 1, 0)
    pub const Y: Self = Self::new(0.0, 1.0, 0.0);

    /// Unit vector along Z axis (0, 0, 1)
    pub const Z: Self = Self::new(0.0, 0.0, 1.0);

    /// One vector (1, 1, 1)
    pub const ONE: Self = Self::new(1.0, 1.0, 1.0);

    /// Creates a new Vec3 from components
    #[inline]
    pub const fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x, y, z }
    }

    /// Creates a Vec3 with all components set to the same value
    #[inline]
    pub const fn splat(v: f32) -> Self {
        Self::new(v, v, v)
    }

    /// Creates a Vec3 from an array
    #[inline]
    pub const fn from_array(arr: [f32; 3]) -> Self {
        Self::new(arr[0], arr[1], arr[2])
    }

    /// Converts the Vec3 to an array
    #[inline]
    pub const fn to_array(self) -> [f32; 3] {
        [self.x, self.y, self.z]
    }

    /// Dot product of two vectors
    #[inline]
    pub fn dot(self, other: Self) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Cross product of two vectors
    #[inline]
    pub fn cross(self, other: Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    /// Squared length of the vector (avoids sqrt)
    #[inline]
    pub fn length_squared(self) -> f32 {
        self.dot(self)
    }

    /// Length (magnitude) of the vector
    #[inline]
    pub fn length(self) -> f32 {
        self.length_squared().sqrt()
    }

    /// Returns a normalized (unit length) version of the vector
    /// Returns zero vector if the input is zero or near-zero
    #[inline]
    pub fn normalize(self) -> Self {
        let len_sq = self.length_squared();
        if len_sq > 1e-10 {
            self / len_sq.sqrt()
        } else {
            Self::ZERO
        }
    }

    /// Returns a normalized vector and its original length
    /// Returns (ZERO, 0.0) if the input is zero or near-zero
    #[inline]
    pub fn normalize_with_length(self) -> (Self, f32) {
        let len_sq = self.length_squared();
        if len_sq > 1e-10 {
            let len = len_sq.sqrt();
            (self / len, len)
        } else {
            (Self::ZERO, 0.0)
        }
    }

    /// Attempts to normalize, returning None if the vector is too small
    #[inline]
    pub fn try_normalize(self) -> Option<Self> {
        let len_sq = self.length_squared();
        if len_sq > 1e-10 {
            Some(self / len_sq.sqrt())
        } else {
            None
        }
    }

    /// Returns true if the vector is approximately zero
    #[inline]
    pub fn is_near_zero(self, epsilon: f32) -> bool {
        self.length_squared() < epsilon * epsilon
    }

    /// Linear interpolation between two vectors
    #[inline]
    pub fn lerp(self, other: Self, t: f32) -> Self {
        self + (other - self) * t
    }

    /// Component-wise minimum
    #[inline]
    pub fn min(self, other: Self) -> Self {
        Self::new(
            self.x.min(other.x),
            self.y.min(other.y),
            self.z.min(other.z),
        )
    }

    /// Component-wise maximum
    #[inline]
    pub fn max(self, other: Self) -> Self {
        Self::new(
            self.x.max(other.x),
            self.y.max(other.y),
            self.z.max(other.z),
        )
    }

    /// Clamps each component to the range [min, max]
    #[inline]
    pub fn clamp(self, min: Self, max: Self) -> Self {
        self.max(min).min(max)
    }

    /// Component-wise absolute value
    #[inline]
    pub fn abs(self) -> Self {
        Self::new(self.x.abs(), self.y.abs(), self.z.abs())
    }

    /// Returns the smallest component
    #[inline]
    pub fn min_element(self) -> f32 {
        self.x.min(self.y).min(self.z)
    }

    /// Returns the largest component
    #[inline]
    pub fn max_element(self) -> f32 {
        self.x.max(self.y).max(self.z)
    }

    /// Component-wise multiplication (Hadamard product)
    #[inline]
    pub fn component_mul(self, other: Self) -> Self {
        Self::new(self.x * other.x, self.y * other.y, self.z * other.z)
    }

    /// Component-wise division
    #[inline]
    pub fn component_div(self, other: Self) -> Self {
        Self::new(self.x / other.x, self.y / other.y, self.z / other.z)
    }

    /// Projects this vector onto another vector
    #[inline]
    pub fn project_onto(self, other: Self) -> Self {
        let other_len_sq = other.length_squared();
        if other_len_sq > 1e-10 {
            other * (self.dot(other) / other_len_sq)
        } else {
            Self::ZERO
        }
    }

    /// Returns the component of this vector perpendicular to another
    #[inline]
    pub fn reject_from(self, other: Self) -> Self {
        self - self.project_onto(other)
    }

    /// Reflects this vector about a normal
    #[inline]
    pub fn reflect(self, normal: Self) -> Self {
        self - normal * (2.0 * self.dot(normal))
    }

    /// Returns the distance between two points
    #[inline]
    pub fn distance(self, other: Self) -> f32 {
        (other - self).length()
    }

    /// Returns the squared distance between two points
    #[inline]
    pub fn distance_squared(self, other: Self) -> f32 {
        (other - self).length_squared()
    }

    /// Returns an arbitrary vector perpendicular to this one
    #[inline]
    pub fn any_perpendicular(self) -> Self {
        let abs_x = self.x.abs();
        let abs_y = self.y.abs();
        let abs_z = self.z.abs();

        // Cross with the axis that is most perpendicular
        let other = if abs_x <= abs_y && abs_x <= abs_z {
            Self::X
        } else if abs_y <= abs_z {
            Self::Y
        } else {
            Self::Z
        };

        self.cross(other).normalize()
    }
}

// Operator overloads

impl Add for Vec3 {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl AddAssign for Vec3 {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl Sub for Vec3 {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl SubAssign for Vec3 {
    #[inline]
    fn sub_assign(&mut self, other: Self) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }
}

impl Mul<f32> for Vec3 {
    type Output = Self;

    #[inline]
    fn mul(self, scalar: f32) -> Self {
        Self::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }
}

impl Mul<Vec3> for f32 {
    type Output = Vec3;

    #[inline]
    fn mul(self, vec: Vec3) -> Vec3 {
        Vec3::new(self * vec.x, self * vec.y, self * vec.z)
    }
}

impl MulAssign<f32> for Vec3 {
    #[inline]
    fn mul_assign(&mut self, scalar: f32) {
        self.x *= scalar;
        self.y *= scalar;
        self.z *= scalar;
    }
}

impl Div<f32> for Vec3 {
    type Output = Self;

    #[inline]
    fn div(self, scalar: f32) -> Self {
        let inv = 1.0 / scalar;
        Self::new(self.x * inv, self.y * inv, self.z * inv)
    }
}

impl DivAssign<f32> for Vec3 {
    #[inline]
    fn div_assign(&mut self, scalar: f32) {
        let inv = 1.0 / scalar;
        self.x *= inv;
        self.y *= inv;
        self.z *= inv;
    }
}

impl Neg for Vec3 {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl Index<usize> for Vec3 {
    type Output = f32;

    #[inline]
    fn index(&self, index: usize) -> &f32 {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Vec3 index out of bounds: {}", index),
        }
    }
}

impl IndexMut<usize> for Vec3 {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut f32 {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("Vec3 index out of bounds: {}", index),
        }
    }
}

impl From<[f32; 3]> for Vec3 {
    #[inline]
    fn from(arr: [f32; 3]) -> Self {
        Self::from_array(arr)
    }
}

impl From<Vec3> for [f32; 3] {
    #[inline]
    fn from(v: Vec3) -> Self {
        v.to_array()
    }
}

impl From<(f32, f32, f32)> for Vec3 {
    #[inline]
    fn from((x, y, z): (f32, f32, f32)) -> Self {
        Self::new(x, y, z)
    }
}

impl From<Vec3> for (f32, f32, f32) {
    #[inline]
    fn from(v: Vec3) -> Self {
        (v.x, v.y, v.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f32 = 1e-6;

    fn approx_eq(a: f32, b: f32) -> bool {
        (a - b).abs() < EPSILON
    }

    fn vec3_approx_eq(a: Vec3, b: Vec3) -> bool {
        approx_eq(a.x, b.x) && approx_eq(a.y, b.y) && approx_eq(a.z, b.z)
    }

    #[test]
    fn test_construction() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        assert_eq!(v.x, 1.0);
        assert_eq!(v.y, 2.0);
        assert_eq!(v.z, 3.0);

        let v2 = Vec3::splat(5.0);
        assert_eq!(v2, Vec3::new(5.0, 5.0, 5.0));
    }

    #[test]
    fn test_dot_product() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        assert!(approx_eq(a.dot(b), 32.0)); // 1*4 + 2*5 + 3*6 = 32
    }

    #[test]
    fn test_cross_product() {
        let x = Vec3::X;
        let y = Vec3::Y;
        let z = x.cross(y);
        assert!(vec3_approx_eq(z, Vec3::Z));

        // Cross product is anti-commutative
        let z2 = y.cross(x);
        assert!(vec3_approx_eq(z2, -Vec3::Z));
    }

    #[test]
    fn test_length() {
        let v = Vec3::new(3.0, 4.0, 0.0);
        assert!(approx_eq(v.length(), 5.0));
        assert!(approx_eq(v.length_squared(), 25.0));
    }

    #[test]
    fn test_normalize() {
        let v = Vec3::new(3.0, 4.0, 0.0);
        let n = v.normalize();
        assert!(approx_eq(n.length(), 1.0));
        assert!(vec3_approx_eq(n, Vec3::new(0.6, 0.8, 0.0)));

        // Zero vector should normalize to zero
        let zero = Vec3::ZERO.normalize();
        assert!(vec3_approx_eq(zero, Vec3::ZERO));
    }

    #[test]
    fn test_lerp() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(10.0, 20.0, 30.0);

        assert!(vec3_approx_eq(a.lerp(b, 0.0), a));
        assert!(vec3_approx_eq(a.lerp(b, 1.0), b));
        assert!(vec3_approx_eq(a.lerp(b, 0.5), Vec3::new(5.0, 10.0, 15.0)));
    }

    #[test]
    fn test_operators() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);

        assert!(vec3_approx_eq(a + b, Vec3::new(5.0, 7.0, 9.0)));
        assert!(vec3_approx_eq(b - a, Vec3::new(3.0, 3.0, 3.0)));
        assert!(vec3_approx_eq(a * 2.0, Vec3::new(2.0, 4.0, 6.0)));
        assert!(vec3_approx_eq(2.0 * a, Vec3::new(2.0, 4.0, 6.0)));
        assert!(vec3_approx_eq(a / 2.0, Vec3::new(0.5, 1.0, 1.5)));
        assert!(vec3_approx_eq(-a, Vec3::new(-1.0, -2.0, -3.0)));
    }

    #[test]
    fn test_project() {
        let v = Vec3::new(3.0, 4.0, 0.0);
        let onto = Vec3::X;
        let projected = v.project_onto(onto);
        assert!(vec3_approx_eq(projected, Vec3::new(3.0, 0.0, 0.0)));
    }

    #[test]
    fn test_reflect() {
        let v = Vec3::new(1.0, -1.0, 0.0);
        let normal = Vec3::Y;
        let reflected = v.reflect(normal);
        assert!(vec3_approx_eq(reflected, Vec3::new(1.0, 1.0, 0.0)));
    }

    #[test]
    fn test_indexing() {
        let mut v = Vec3::new(1.0, 2.0, 3.0);
        assert_eq!(v[0], 1.0);
        assert_eq!(v[1], 2.0);
        assert_eq!(v[2], 3.0);

        v[1] = 5.0;
        assert_eq!(v[1], 5.0);
    }
}
