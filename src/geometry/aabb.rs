use crate::math::Vec3;

/// An axis-aligned bounding box defined by minimum and maximum points.
///
/// Used for broad-phase collision detection and spatial queries.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb {
    /// Minimum corner (smallest x, y, z values)
    pub min: Vec3,
    /// Maximum corner (largest x, y, z values)
    pub max: Vec3,
}

impl Default for Aabb {
    fn default() -> Self {
        Self::EMPTY
    }
}

impl Aabb {
    /// An empty AABB that contains no points
    pub const EMPTY: Self = Self {
        min: Vec3::new(f32::INFINITY, f32::INFINITY, f32::INFINITY),
        max: Vec3::new(f32::NEG_INFINITY, f32::NEG_INFINITY, f32::NEG_INFINITY),
    };

    /// An infinitely large AABB that contains all points
    pub const INFINITE: Self = Self {
        min: Vec3::new(f32::NEG_INFINITY, f32::NEG_INFINITY, f32::NEG_INFINITY),
        max: Vec3::new(f32::INFINITY, f32::INFINITY, f32::INFINITY),
    };

    /// Creates an AABB from minimum and maximum points
    #[inline]
    pub const fn new(min: Vec3, max: Vec3) -> Self {
        Self { min, max }
    }

    /// Creates an AABB from center and half-extents
    #[inline]
    pub fn from_center_half_extents(center: Vec3, half_extents: Vec3) -> Self {
        Self {
            min: center - half_extents,
            max: center + half_extents,
        }
    }

    /// Creates an AABB that contains a single point
    #[inline]
    pub fn from_point(point: Vec3) -> Self {
        Self {
            min: point,
            max: point,
        }
    }

    /// Creates an AABB from an array of points
    #[inline]
    pub fn from_points(points: &[Vec3]) -> Self {
        let mut aabb = Self::EMPTY;
        for &point in points {
            aabb = aabb.expand_to_include(point);
        }
        aabb
    }

    /// Returns the center of the AABB
    #[inline]
    pub fn center(self) -> Vec3 {
        (self.min + self.max) * 0.5
    }

    /// Returns the half-extents (half the size in each dimension)
    #[inline]
    pub fn half_extents(self) -> Vec3 {
        (self.max - self.min) * 0.5
    }

    /// Returns the full size (extents) of the AABB
    #[inline]
    pub fn size(self) -> Vec3 {
        self.max - self.min
    }

    /// Returns the volume of the AABB
    #[inline]
    pub fn volume(self) -> f32 {
        let size = self.size();
        size.x * size.y * size.z
    }

    /// Returns the surface area of the AABB
    #[inline]
    pub fn surface_area(self) -> f32 {
        let size = self.size();
        2.0 * (size.x * size.y + size.y * size.z + size.z * size.x)
    }

    /// Returns true if this AABB is valid (min <= max in all dimensions)
    #[inline]
    pub fn is_valid(self) -> bool {
        self.min.x <= self.max.x && self.min.y <= self.max.y && self.min.z <= self.max.z
    }

    /// Returns true if this AABB is empty
    #[inline]
    pub fn is_empty(self) -> bool {
        self.min.x > self.max.x || self.min.y > self.max.y || self.min.z > self.max.z
    }

    /// Returns true if this AABB contains the given point
    #[inline]
    pub fn contains_point(self, point: Vec3) -> bool {
        point.x >= self.min.x
            && point.x <= self.max.x
            && point.y >= self.min.y
            && point.y <= self.max.y
            && point.z >= self.min.z
            && point.z <= self.max.z
    }

    /// Returns true if this AABB fully contains another AABB
    #[inline]
    pub fn contains_aabb(self, other: Self) -> bool {
        self.min.x <= other.min.x
            && self.max.x >= other.max.x
            && self.min.y <= other.min.y
            && self.max.y >= other.max.y
            && self.min.z <= other.min.z
            && self.max.z >= other.max.z
    }

    /// Returns true if this AABB intersects another AABB
    #[inline]
    pub fn intersects(self, other: Self) -> bool {
        self.min.x <= other.max.x
            && self.max.x >= other.min.x
            && self.min.y <= other.max.y
            && self.max.y >= other.min.y
            && self.min.z <= other.max.z
            && self.max.z >= other.min.z
    }

    /// Returns the intersection of two AABBs, or EMPTY if they don't intersect
    #[inline]
    pub fn intersection(self, other: Self) -> Self {
        let min = self.min.max(other.min);
        let max = self.max.min(other.max);
        Self { min, max }
    }

    /// Returns a new AABB that is the union of this and another AABB
    #[inline]
    pub fn union(self, other: Self) -> Self {
        Self {
            min: self.min.min(other.min),
            max: self.max.max(other.max),
        }
    }

    /// Returns a new AABB expanded to include a point
    #[inline]
    pub fn expand_to_include(self, point: Vec3) -> Self {
        Self {
            min: self.min.min(point),
            max: self.max.max(point),
        }
    }

    /// Returns a new AABB expanded by a margin in all directions
    #[inline]
    pub fn expand(self, margin: f32) -> Self {
        let m = Vec3::splat(margin);
        Self {
            min: self.min - m,
            max: self.max + m,
        }
    }

    /// Returns the closest point on the AABB to the given point
    #[inline]
    pub fn closest_point(self, point: Vec3) -> Vec3 {
        point.clamp(self.min, self.max)
    }

    /// Returns the squared distance from a point to the AABB
    #[inline]
    pub fn distance_squared_to_point(self, point: Vec3) -> f32 {
        let closest = self.closest_point(point);
        closest.distance_squared(point)
    }

    /// Returns the distance from a point to the AABB
    #[inline]
    pub fn distance_to_point(self, point: Vec3) -> f32 {
        self.distance_squared_to_point(point).sqrt()
    }

    /// Tests intersection with a ray.
    /// Returns Some((t_min, t_max)) if the ray intersects, where t_min and t_max
    /// are the entry and exit distances along the ray.
    #[inline]
    pub fn ray_intersection(self, origin: Vec3, direction: Vec3) -> Option<(f32, f32)> {
        let inv_dir = Vec3::new(1.0 / direction.x, 1.0 / direction.y, 1.0 / direction.z);

        let t1 = (self.min.x - origin.x) * inv_dir.x;
        let t2 = (self.max.x - origin.x) * inv_dir.x;
        let t3 = (self.min.y - origin.y) * inv_dir.y;
        let t4 = (self.max.y - origin.y) * inv_dir.y;
        let t5 = (self.min.z - origin.z) * inv_dir.z;
        let t6 = (self.max.z - origin.z) * inv_dir.z;

        let t_min = t1.min(t2).max(t3.min(t4)).max(t5.min(t6));
        let t_max = t1.max(t2).min(t3.max(t4)).min(t5.max(t6));

        if t_max >= t_min && t_max >= 0.0 {
            Some((t_min.max(0.0), t_max))
        } else {
            None
        }
    }

    /// Returns the 8 corners of the AABB
    #[inline]
    pub fn corners(self) -> [Vec3; 8] {
        [
            Vec3::new(self.min.x, self.min.y, self.min.z),
            Vec3::new(self.max.x, self.min.y, self.min.z),
            Vec3::new(self.min.x, self.max.y, self.min.z),
            Vec3::new(self.max.x, self.max.y, self.min.z),
            Vec3::new(self.min.x, self.min.y, self.max.z),
            Vec3::new(self.max.x, self.min.y, self.max.z),
            Vec3::new(self.min.x, self.max.y, self.max.z),
            Vec3::new(self.max.x, self.max.y, self.max.z),
        ]
    }

    /// Returns the longest axis (0 = x, 1 = y, 2 = z)
    #[inline]
    pub fn longest_axis(self) -> usize {
        let size = self.size();
        if size.x >= size.y && size.x >= size.z {
            0
        } else if size.y >= size.z {
            1
        } else {
            2
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_creation() {
        let aabb = Aabb::new(Vec3::new(-1.0, -2.0, -3.0), Vec3::new(1.0, 2.0, 3.0));
        assert_eq!(aabb.center(), Vec3::ZERO);
        assert_eq!(aabb.half_extents(), Vec3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn test_from_center_half_extents() {
        let aabb = Aabb::from_center_half_extents(Vec3::new(1.0, 2.0, 3.0), Vec3::new(1.0, 1.0, 1.0));
        assert_eq!(aabb.min, Vec3::new(0.0, 1.0, 2.0));
        assert_eq!(aabb.max, Vec3::new(2.0, 3.0, 4.0));
    }

    #[test]
    fn test_contains_point() {
        let aabb = Aabb::new(Vec3::ZERO, Vec3::ONE);
        assert!(aabb.contains_point(Vec3::new(0.5, 0.5, 0.5)));
        assert!(aabb.contains_point(Vec3::ZERO));
        assert!(aabb.contains_point(Vec3::ONE));
        assert!(!aabb.contains_point(Vec3::new(2.0, 0.5, 0.5)));
    }

    #[test]
    fn test_intersects() {
        let a = Aabb::new(Vec3::ZERO, Vec3::ONE);
        let b = Aabb::new(Vec3::new(0.5, 0.5, 0.5), Vec3::new(1.5, 1.5, 1.5));
        let c = Aabb::new(Vec3::new(2.0, 0.0, 0.0), Vec3::new(3.0, 1.0, 1.0));

        assert!(a.intersects(b));
        assert!(b.intersects(a));
        assert!(!a.intersects(c));
    }

    #[test]
    fn test_union() {
        let a = Aabb::new(Vec3::ZERO, Vec3::ONE);
        let b = Aabb::new(Vec3::new(2.0, 2.0, 2.0), Vec3::new(3.0, 3.0, 3.0));
        let u = a.union(b);

        assert_eq!(u.min, Vec3::ZERO);
        assert_eq!(u.max, Vec3::new(3.0, 3.0, 3.0));
    }

    #[test]
    fn test_expand() {
        let aabb = Aabb::new(Vec3::ZERO, Vec3::ONE);
        let expanded = aabb.expand(0.5);

        assert_eq!(expanded.min, Vec3::new(-0.5, -0.5, -0.5));
        assert_eq!(expanded.max, Vec3::new(1.5, 1.5, 1.5));
    }

    #[test]
    fn test_volume_and_surface_area() {
        let aabb = Aabb::new(Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0));
        assert_eq!(aabb.volume(), 24.0);
        assert_eq!(aabb.surface_area(), 52.0); // 2*(2*3 + 3*4 + 4*2) = 2*(6+12+8) = 52
    }

    #[test]
    fn test_closest_point() {
        let aabb = Aabb::new(Vec3::ZERO, Vec3::ONE);

        // Point inside
        assert_eq!(
            aabb.closest_point(Vec3::new(0.5, 0.5, 0.5)),
            Vec3::new(0.5, 0.5, 0.5)
        );

        // Point outside
        assert_eq!(
            aabb.closest_point(Vec3::new(2.0, 0.5, 0.5)),
            Vec3::new(1.0, 0.5, 0.5)
        );
    }

    #[test]
    fn test_ray_intersection() {
        let aabb = Aabb::new(Vec3::ZERO, Vec3::ONE);

        // Ray that hits
        let result = aabb.ray_intersection(Vec3::new(-1.0, 0.5, 0.5), Vec3::X);
        assert!(result.is_some());
        let (t_min, t_max) = result.unwrap();
        assert!((t_min - 1.0).abs() < 1e-6);
        assert!((t_max - 2.0).abs() < 1e-6);

        // Ray that misses
        let result = aabb.ray_intersection(Vec3::new(-1.0, 2.0, 0.5), Vec3::X);
        assert!(result.is_none());
    }

    #[test]
    fn test_from_points() {
        let points = vec![
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(-1.0, 0.0, 5.0),
            Vec3::new(0.0, -2.0, 1.0),
        ];
        let aabb = Aabb::from_points(&points);

        assert_eq!(aabb.min, Vec3::new(-1.0, -2.0, 1.0));
        assert_eq!(aabb.max, Vec3::new(1.0, 2.0, 5.0));
    }
}
