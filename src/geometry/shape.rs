use crate::math::{Mat3, Transform, Vec3};

use super::aabb::Aabb;

/// The type of collision shape
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShapeType {
    Sphere,
    Box,
    Capsule,
}

/// A collision shape that can be attached to rigid bodies.
#[derive(Debug, Clone, Copy)]
pub enum Shape {
    /// A sphere defined by its radius
    Sphere(Sphere),
    /// A box (cuboid) defined by half-extents
    Box(BoxShape),
    /// A capsule defined by radius and half-height
    Capsule(Capsule),
}

impl Shape {
    /// Creates a sphere shape
    #[inline]
    pub fn sphere(radius: f32) -> Self {
        Self::Sphere(Sphere::new(radius))
    }

    /// Creates a box shape from half-extents
    #[inline]
    pub fn cuboid(half_extents: Vec3) -> Self {
        Self::Box(BoxShape::new(half_extents))
    }

    /// Creates a box shape from full dimensions
    #[inline]
    pub fn cuboid_from_size(size: Vec3) -> Self {
        Self::Box(BoxShape::from_size(size))
    }

    /// Creates a capsule shape
    #[inline]
    pub fn capsule(radius: f32, half_height: f32) -> Self {
        Self::Capsule(Capsule::new(radius, half_height))
    }

    /// Returns the shape type
    #[inline]
    pub fn shape_type(&self) -> ShapeType {
        match self {
            Shape::Sphere(_) => ShapeType::Sphere,
            Shape::Box(_) => ShapeType::Box,
            Shape::Capsule(_) => ShapeType::Capsule,
        }
    }

    /// Computes the AABB of this shape in local space
    #[inline]
    pub fn local_aabb(&self) -> Aabb {
        match self {
            Shape::Sphere(s) => s.aabb(),
            Shape::Box(b) => b.aabb(),
            Shape::Capsule(c) => c.aabb(),
        }
    }

    /// Computes the AABB of this shape given a world transform
    #[inline]
    pub fn world_aabb(&self, transform: Transform) -> Aabb {
        match self {
            Shape::Sphere(s) => s.world_aabb(transform),
            Shape::Box(b) => b.world_aabb(transform),
            Shape::Capsule(c) => c.world_aabb(transform),
        }
    }

    /// Returns the support point in the given direction (local space)
    /// The support point is the point on the shape furthest in the given direction
    #[inline]
    pub fn support(&self, direction: Vec3) -> Vec3 {
        match self {
            Shape::Sphere(s) => s.support(direction),
            Shape::Box(b) => b.support(direction),
            Shape::Capsule(c) => c.support(direction),
        }
    }

    /// Returns the support point in world space given a transform
    #[inline]
    pub fn support_world(&self, transform: Transform, direction: Vec3) -> Vec3 {
        // Transform direction to local space
        let local_dir = transform.inverse_transform_vector(direction);
        // Get support in local space
        let local_support = self.support(local_dir);
        // Transform back to world space
        transform.transform_point(local_support)
    }

    /// Computes mass properties (mass, center of mass, inertia tensor) given density
    #[inline]
    pub fn mass_properties(&self, density: f32) -> MassProperties {
        match self {
            Shape::Sphere(s) => s.mass_properties(density),
            Shape::Box(b) => b.mass_properties(density),
            Shape::Capsule(c) => c.mass_properties(density),
        }
    }
}

/// Mass properties of a shape
#[derive(Debug, Clone, Copy)]
pub struct MassProperties {
    /// Total mass
    pub mass: f32,
    /// Center of mass in local coordinates (usually zero for symmetric shapes)
    pub center_of_mass: Vec3,
    /// Inertia tensor in local coordinates (diagonal for primitive shapes)
    pub inertia: Mat3,
}

impl MassProperties {
    /// Creates mass properties with zero mass (for static bodies)
    pub const ZERO: Self = Self {
        mass: 0.0,
        center_of_mass: Vec3::ZERO,
        inertia: Mat3::ZERO,
    };

    /// Creates mass properties from mass and inertia diagonal
    #[inline]
    pub fn new(mass: f32, inertia_diagonal: Vec3) -> Self {
        Self {
            mass,
            center_of_mass: Vec3::ZERO,
            inertia: Mat3::from_diagonal(inertia_diagonal),
        }
    }

    /// Returns the inverse mass (0 for infinite mass)
    #[inline]
    pub fn inv_mass(&self) -> f32 {
        if self.mass > 0.0 {
            1.0 / self.mass
        } else {
            0.0
        }
    }

    /// Returns the inverse inertia tensor
    #[inline]
    pub fn inv_inertia(&self) -> Mat3 {
        self.inertia.try_inverse().unwrap_or(Mat3::ZERO)
    }
}

/// A sphere collision shape
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere {
    pub radius: f32,
}

impl Sphere {
    /// Creates a new sphere with the given radius
    #[inline]
    pub fn new(radius: f32) -> Self {
        Self { radius }
    }

    /// Returns the AABB of this sphere in local space
    #[inline]
    pub fn aabb(&self) -> Aabb {
        let r = Vec3::splat(self.radius);
        Aabb::new(-r, r)
    }

    /// Returns the AABB of this sphere given a world transform
    #[inline]
    pub fn world_aabb(&self, transform: Transform) -> Aabb {
        let r = Vec3::splat(self.radius);
        Aabb::new(transform.position - r, transform.position + r)
    }

    /// Returns the support point in the given direction
    #[inline]
    pub fn support(&self, direction: Vec3) -> Vec3 {
        direction.normalize() * self.radius
    }

    /// Computes mass properties given density
    #[inline]
    pub fn mass_properties(&self, density: f32) -> MassProperties {
        let r = self.radius;
        let volume = (4.0 / 3.0) * std::f32::consts::PI * r * r * r;
        let mass = volume * density;
        let i = (2.0 / 5.0) * mass * r * r;
        MassProperties::new(mass, Vec3::splat(i))
    }

    /// Returns the volume of the sphere
    #[inline]
    pub fn volume(&self) -> f32 {
        (4.0 / 3.0) * std::f32::consts::PI * self.radius * self.radius * self.radius
    }
}

/// A box (cuboid) collision shape
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoxShape {
    /// Half-extents (half the size in each dimension)
    pub half_extents: Vec3,
}

impl BoxShape {
    /// Creates a new box with the given half-extents
    #[inline]
    pub fn new(half_extents: Vec3) -> Self {
        Self { half_extents }
    }

    /// Creates a box from full size dimensions
    #[inline]
    pub fn from_size(size: Vec3) -> Self {
        Self {
            half_extents: size * 0.5,
        }
    }

    /// Returns the full size of the box
    #[inline]
    pub fn size(&self) -> Vec3 {
        self.half_extents * 2.0
    }

    /// Returns the AABB of this box in local space
    #[inline]
    pub fn aabb(&self) -> Aabb {
        Aabb::new(-self.half_extents, self.half_extents)
    }

    /// Returns the AABB of this box given a world transform
    #[inline]
    pub fn world_aabb(&self, transform: Transform) -> Aabb {
        let rot = transform.rotation_matrix();

        // Compute the world-space extent by summing absolute values of rotated axes
        let abs_rot = Mat3::from_cols(rot.col(0).abs(), rot.col(1).abs(), rot.col(2).abs());
        let world_half_extents = abs_rot * self.half_extents;

        Aabb::new(
            transform.position - world_half_extents,
            transform.position + world_half_extents,
        )
    }

    /// Returns the support point in the given direction
    #[inline]
    pub fn support(&self, direction: Vec3) -> Vec3 {
        Vec3::new(
            if direction.x >= 0.0 {
                self.half_extents.x
            } else {
                -self.half_extents.x
            },
            if direction.y >= 0.0 {
                self.half_extents.y
            } else {
                -self.half_extents.y
            },
            if direction.z >= 0.0 {
                self.half_extents.z
            } else {
                -self.half_extents.z
            },
        )
    }

    /// Computes mass properties given density
    #[inline]
    pub fn mass_properties(&self, density: f32) -> MassProperties {
        let size = self.size();
        let volume = size.x * size.y * size.z;
        let mass = volume * density;

        let x2 = size.x * size.x;
        let y2 = size.y * size.y;
        let z2 = size.z * size.z;

        let i = Vec3::new(
            (1.0 / 12.0) * mass * (y2 + z2),
            (1.0 / 12.0) * mass * (x2 + z2),
            (1.0 / 12.0) * mass * (x2 + y2),
        );

        MassProperties::new(mass, i)
    }

    /// Returns the volume of the box
    #[inline]
    pub fn volume(&self) -> f32 {
        8.0 * self.half_extents.x * self.half_extents.y * self.half_extents.z
    }

    /// Returns the 8 vertices of the box in local space
    #[inline]
    pub fn vertices(&self) -> [Vec3; 8] {
        let h = self.half_extents;
        [
            Vec3::new(-h.x, -h.y, -h.z),
            Vec3::new(h.x, -h.y, -h.z),
            Vec3::new(-h.x, h.y, -h.z),
            Vec3::new(h.x, h.y, -h.z),
            Vec3::new(-h.x, -h.y, h.z),
            Vec3::new(h.x, -h.y, h.z),
            Vec3::new(-h.x, h.y, h.z),
            Vec3::new(h.x, h.y, h.z),
        ]
    }
}

/// A capsule collision shape (cylinder with hemispherical caps)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Capsule {
    /// Radius of the cylinder and hemispherical caps
    pub radius: f32,
    /// Half-height of the cylindrical segment (total height = 2 * half_height + 2 * radius)
    pub half_height: f32,
}

impl Capsule {
    /// Creates a new capsule with the given radius and half-height
    /// The capsule is oriented along the Y axis
    #[inline]
    pub fn new(radius: f32, half_height: f32) -> Self {
        Self { radius, half_height }
    }

    /// Returns the total height of the capsule (including caps)
    #[inline]
    pub fn total_height(&self) -> f32 {
        2.0 * self.half_height + 2.0 * self.radius
    }

    /// Returns the AABB of this capsule in local space
    #[inline]
    pub fn aabb(&self) -> Aabb {
        let r = self.radius;
        let h = self.half_height + r;
        Aabb::new(Vec3::new(-r, -h, -r), Vec3::new(r, h, r))
    }

    /// Returns the AABB of this capsule given a world transform
    #[inline]
    pub fn world_aabb(&self, transform: Transform) -> Aabb {
        // The capsule's axis is along local Y
        let local_axis = Vec3::Y;
        let world_axis = transform.transform_vector(local_axis);

        // End points of the center line
        let p0 = transform.position - world_axis * self.half_height;
        let p1 = transform.position + world_axis * self.half_height;

        // Create AABB around the two sphere centers and expand by radius
        let min = p0.min(p1) - Vec3::splat(self.radius);
        let max = p0.max(p1) + Vec3::splat(self.radius);

        Aabb::new(min, max)
    }

    /// Returns the support point in the given direction
    #[inline]
    pub fn support(&self, direction: Vec3) -> Vec3 {
        // Project direction onto capsule axis (Y)
        let center = if direction.y >= 0.0 {
            Vec3::new(0.0, self.half_height, 0.0)
        } else {
            Vec3::new(0.0, -self.half_height, 0.0)
        };

        // Add sphere support at that center
        center + direction.normalize() * self.radius
    }

    /// Computes mass properties given density
    #[inline]
    pub fn mass_properties(&self, density: f32) -> MassProperties {
        let r = self.radius;
        let h = self.half_height * 2.0; // Full cylinder height

        // Cylinder volume + sphere volume
        let cylinder_volume = std::f32::consts::PI * r * r * h;
        let sphere_volume = (4.0 / 3.0) * std::f32::consts::PI * r * r * r;
        let volume = cylinder_volume + sphere_volume;
        let mass = volume * density;

        // Approximate inertia (treating as cylinder with hemispherical caps)
        let cylinder_mass = cylinder_volume * density;
        let sphere_mass = sphere_volume * density;

        // Cylinder inertia (aligned along Y)
        let i_cyl_y = 0.5 * cylinder_mass * r * r;
        let i_cyl_xz = (1.0 / 12.0) * cylinder_mass * (3.0 * r * r + h * h);

        // Sphere inertia (at distance half_height from center)
        let i_sphere = (2.0 / 5.0) * sphere_mass * r * r;
        // Parallel axis theorem for the two hemispheres
        let i_sphere_offset = i_sphere + sphere_mass * self.half_height * self.half_height;

        let i_y = i_cyl_y + i_sphere;
        let i_xz = i_cyl_xz + i_sphere_offset;

        MassProperties::new(mass, Vec3::new(i_xz, i_y, i_xz))
    }

    /// Returns the volume of the capsule
    #[inline]
    pub fn volume(&self) -> f32 {
        let r = self.radius;
        let h = self.half_height * 2.0;
        std::f32::consts::PI * r * r * h + (4.0 / 3.0) * std::f32::consts::PI * r * r * r
    }

    /// Returns the two endpoints of the capsule's center line
    #[inline]
    pub fn segment(&self) -> (Vec3, Vec3) {
        (
            Vec3::new(0.0, -self.half_height, 0.0),
            Vec3::new(0.0, self.half_height, 0.0),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Quat;
    use std::f32::consts::PI;

    const EPSILON: f32 = 1e-5;

    fn approx_eq(a: f32, b: f32) -> bool {
        (a - b).abs() < EPSILON
    }

    #[test]
    fn test_sphere_support() {
        let sphere = Sphere::new(2.0);

        let support = sphere.support(Vec3::X);
        assert!(approx_eq(support.x, 2.0));
        assert!(approx_eq(support.y, 0.0));
        assert!(approx_eq(support.z, 0.0));

        let support = sphere.support(Vec3::new(1.0, 1.0, 0.0));
        let expected_len = 2.0;
        assert!(approx_eq(support.length(), expected_len));
    }

    #[test]
    fn test_box_support() {
        let b = BoxShape::new(Vec3::new(1.0, 2.0, 3.0));

        assert_eq!(b.support(Vec3::X), Vec3::new(1.0, 2.0, 3.0));
        assert_eq!(b.support(-Vec3::X), Vec3::new(-1.0, 2.0, 3.0));
        assert_eq!(b.support(-Vec3::ONE), Vec3::new(-1.0, -2.0, -3.0));
    }

    #[test]
    fn test_capsule_support() {
        let capsule = Capsule::new(1.0, 2.0);

        // Support along Y should be at top cap
        let support = capsule.support(Vec3::Y);
        assert!(approx_eq(support.y, 3.0)); // half_height + radius

        // Support along -Y should be at bottom cap
        let support = capsule.support(-Vec3::Y);
        assert!(approx_eq(support.y, -3.0));
    }

    #[test]
    fn test_sphere_mass() {
        let sphere = Sphere::new(1.0);
        let props = sphere.mass_properties(1.0);

        let expected_volume = (4.0 / 3.0) * PI;
        assert!(approx_eq(props.mass, expected_volume));
    }

    #[test]
    fn test_box_mass() {
        let b = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let props = b.mass_properties(1.0);

        // Volume = 8 (cube with side 2)
        assert!(approx_eq(props.mass, 8.0));
    }

    #[test]
    fn test_box_world_aabb() {
        let b = BoxShape::new(Vec3::new(1.0, 2.0, 3.0));

        // No rotation - AABB should match local AABB translated
        let t = Transform::from_position(Vec3::new(10.0, 0.0, 0.0));
        let aabb = b.world_aabb(t);
        assert_eq!(aabb.min, Vec3::new(9.0, -2.0, -3.0));
        assert_eq!(aabb.max, Vec3::new(11.0, 2.0, 3.0));

        // 45 degree rotation around Z
        let t = Transform::from_rotation(Quat::from_axis_angle(Vec3::Z, PI / 4.0));
        let aabb = b.world_aabb(t);
        // The rotated box should have a larger AABB
        let local_aabb = b.aabb();
        assert!(aabb.size().x > local_aabb.size().x || aabb.size().y > local_aabb.size().y);
    }

    #[test]
    fn test_shape_enum() {
        let sphere = Shape::sphere(2.0);
        assert_eq!(sphere.shape_type(), ShapeType::Sphere);

        let box_shape = Shape::cuboid(Vec3::ONE);
        assert_eq!(box_shape.shape_type(), ShapeType::Box);

        let capsule = Shape::capsule(1.0, 2.0);
        assert_eq!(capsule.shape_type(), ShapeType::Capsule);
    }

    #[test]
    fn test_support_world() {
        let shape = Shape::sphere(1.0);
        let transform = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));

        let support = shape.support_world(transform, Vec3::X);
        assert!(approx_eq(support.x, 6.0)); // 5 + 1
    }
}
