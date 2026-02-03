mod mat3;
mod quat;
mod transform;
mod vec3;

pub use mat3::Mat3;
pub use quat::Quat;
pub use transform::{Isometry, Transform};
pub use vec3::Vec3;

/// Common math constants and utilities
pub mod consts {
    /// A small epsilon value for floating point comparisons
    pub const EPSILON: f32 = 1e-6;

    /// Pi constant
    pub const PI: f32 = std::f32::consts::PI;

    /// Two times Pi
    pub const TAU: f32 = std::f32::consts::TAU;

    /// Pi divided by 2
    pub const FRAC_PI_2: f32 = std::f32::consts::FRAC_PI_2;

    /// Pi divided by 4
    pub const FRAC_PI_4: f32 = std::f32::consts::FRAC_PI_4;
}

/// Utility functions
pub mod utils {
    /// Clamps a value to the range [min, max]
    #[inline]
    pub fn clamp(value: f32, min: f32, max: f32) -> f32 {
        value.max(min).min(max)
    }

    /// Linear interpolation between two values
    #[inline]
    pub fn lerp(a: f32, b: f32, t: f32) -> f32 {
        a + (b - a) * t
    }

    /// Returns true if two floats are approximately equal
    #[inline]
    pub fn approx_eq(a: f32, b: f32, epsilon: f32) -> bool {
        (a - b).abs() < epsilon
    }

    /// Converts degrees to radians
    #[inline]
    pub fn deg_to_rad(degrees: f32) -> f32 {
        degrees * (std::f32::consts::PI / 180.0)
    }

    /// Converts radians to degrees
    #[inline]
    pub fn rad_to_deg(radians: f32) -> f32 {
        radians * (180.0 / std::f32::consts::PI)
    }
}
