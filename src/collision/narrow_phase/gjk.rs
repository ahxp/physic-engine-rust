use crate::geometry::Shape;
use crate::math::{Transform, Vec3};

/// Maximum iterations for GJK algorithm
const GJK_MAX_ITERATIONS: usize = 64;

/// Tolerance for GJK convergence
const GJK_TOLERANCE: f32 = 1e-6;

/// Result of a GJK query
#[derive(Debug, Clone)]
pub enum GjkResult {
    /// Shapes are intersecting, simplex can be used for EPA
    Intersecting(Simplex),
    /// Shapes are not intersecting, includes closest points and distance
    Separated {
        closest_a: Vec3,
        closest_b: Vec3,
        distance: f32,
    },
}

impl GjkResult {
    /// Returns true if shapes are intersecting
    pub fn is_intersecting(&self) -> bool {
        matches!(self, GjkResult::Intersecting(_))
    }
}

/// A simplex used in the GJK algorithm
/// Stores up to 4 points (vertices) and their support points on each shape
#[derive(Debug, Clone)]
pub struct Simplex {
    /// Points in the Minkowski difference (A - B space)
    pub points: [Vec3; 4],
    /// Support points on shape A
    pub support_a: [Vec3; 4],
    /// Support points on shape B
    pub support_b: [Vec3; 4],
    /// Number of points in the simplex (1-4)
    pub size: usize,
}

impl Default for Simplex {
    fn default() -> Self {
        Self::new()
    }
}

impl Simplex {
    /// Creates an empty simplex
    pub fn new() -> Self {
        Self {
            points: [Vec3::ZERO; 4],
            support_a: [Vec3::ZERO; 4],
            support_b: [Vec3::ZERO; 4],
            size: 0,
        }
    }

    /// Adds a point to the simplex
    pub fn push(&mut self, point: Vec3, support_a: Vec3, support_b: Vec3) {
        debug_assert!(self.size < 4);
        self.points[self.size] = point;
        self.support_a[self.size] = support_a;
        self.support_b[self.size] = support_b;
        self.size += 1;
    }

    /// Removes point at index and shifts remaining points
    pub fn remove(&mut self, index: usize) {
        debug_assert!(index < self.size);
        for i in index..self.size - 1 {
            self.points[i] = self.points[i + 1];
            self.support_a[i] = self.support_a[i + 1];
            self.support_b[i] = self.support_b[i + 1];
        }
        self.size -= 1;
    }

    /// Gets the last point added
    pub fn last(&self) -> Vec3 {
        debug_assert!(self.size > 0);
        self.points[self.size - 1]
    }

    /// Gets point at index
    pub fn get(&self, index: usize) -> Vec3 {
        debug_assert!(index < self.size);
        self.points[index]
    }
}

/// Computes the support point in the Minkowski difference (A - B)
fn support(
    shape_a: &Shape,
    transform_a: Transform,
    shape_b: &Shape,
    transform_b: Transform,
    direction: Vec3,
) -> (Vec3, Vec3, Vec3) {
    let support_a = shape_a.support_world(transform_a, direction);
    let support_b = shape_b.support_world(transform_b, -direction);
    let point = support_a - support_b;
    (point, support_a, support_b)
}

/// Performs the GJK algorithm to detect collision between two shapes
pub fn gjk(
    shape_a: &Shape,
    transform_a: Transform,
    shape_b: &Shape,
    transform_b: Transform,
) -> GjkResult {
    let mut simplex = Simplex::new();

    // Initial direction from A to B
    let mut direction = transform_b.position - transform_a.position;
    if direction.is_near_zero(GJK_TOLERANCE) {
        direction = Vec3::X;
    }

    // First support point
    let (point, sup_a, sup_b) = support(shape_a, transform_a, shape_b, transform_b, direction);
    simplex.push(point, sup_a, sup_b);

    // New direction towards origin
    direction = -point;

    for _ in 0..GJK_MAX_ITERATIONS {
        if direction.is_near_zero(GJK_TOLERANCE) {
            // Origin is on the simplex boundary
            return GjkResult::Intersecting(simplex);
        }

        let (point, sup_a, sup_b) = support(shape_a, transform_a, shape_b, transform_b, direction);

        // Check if we passed the origin
        if point.dot(direction) < GJK_TOLERANCE {
            // Did not pass origin, shapes are separated
            return compute_closest_points(&simplex, direction);
        }

        simplex.push(point, sup_a, sup_b);

        // Update simplex and direction
        if do_simplex(&mut simplex, &mut direction) {
            return GjkResult::Intersecting(simplex);
        }
    }

    // Max iterations reached, assume separated
    compute_closest_points(&simplex, direction)
}

/// Updates the simplex and direction, returns true if origin is enclosed
fn do_simplex(simplex: &mut Simplex, direction: &mut Vec3) -> bool {
    match simplex.size {
        2 => line_case(simplex, direction),
        3 => triangle_case(simplex, direction),
        4 => tetrahedron_case(simplex, direction),
        _ => false,
    }
}

/// Handles the line case (2 points)
fn line_case(simplex: &mut Simplex, direction: &mut Vec3) -> bool {
    let a = simplex.get(1); // Most recently added
    let b = simplex.get(0);

    let ab = b - a;
    let ao = -a;

    if ab.dot(ao) > 0.0 {
        // Origin is between A and B
        *direction = ab.cross(ao).cross(ab);
        if direction.is_near_zero(GJK_TOLERANCE) {
            // Origin is on the line segment
            *direction = ab.any_perpendicular();
        }
    } else {
        // Origin is beyond A
        simplex.remove(0);
        *direction = ao;
    }

    false
}

/// Handles the triangle case (3 points)
fn triangle_case(simplex: &mut Simplex, direction: &mut Vec3) -> bool {
    let a = simplex.get(2); // Most recently added
    let b = simplex.get(1);
    let c = simplex.get(0);

    let ab = b - a;
    let ac = c - a;
    let ao = -a;

    let abc = ab.cross(ac); // Triangle normal

    // Check if origin is outside edge AB
    let ab_perp = ab.cross(abc);
    if ab_perp.dot(ao) > 0.0 {
        // Origin is outside AB
        simplex.remove(0); // Remove C
        *direction = ab.cross(ao).cross(ab);
        if direction.is_near_zero(GJK_TOLERANCE) {
            *direction = ab.any_perpendicular();
        }
        return false;
    }

    // Check if origin is outside edge AC
    let ac_perp = abc.cross(ac);
    if ac_perp.dot(ao) > 0.0 {
        // Origin is outside AC
        simplex.remove(1); // Remove B
        *direction = ac.cross(ao).cross(ac);
        if direction.is_near_zero(GJK_TOLERANCE) {
            *direction = ac.any_perpendicular();
        }
        return false;
    }

    // Origin is inside the triangle prism
    if abc.dot(ao) > 0.0 {
        // Origin is above triangle
        *direction = abc;
    } else {
        // Origin is below triangle, flip winding
        let temp_point = simplex.points[0];
        let temp_a = simplex.support_a[0];
        let temp_b = simplex.support_b[0];
        simplex.points[0] = simplex.points[1];
        simplex.support_a[0] = simplex.support_a[1];
        simplex.support_b[0] = simplex.support_b[1];
        simplex.points[1] = temp_point;
        simplex.support_a[1] = temp_a;
        simplex.support_b[1] = temp_b;
        *direction = -abc;
    }

    false
}

/// Handles the tetrahedron case (4 points)
fn tetrahedron_case(simplex: &mut Simplex, direction: &mut Vec3) -> bool {
    let a = simplex.get(3); // Most recently added
    let b = simplex.get(2);
    let c = simplex.get(1);
    let d = simplex.get(0);

    let ab = b - a;
    let ac = c - a;
    let ad = d - a;
    let ao = -a;

    // Check each face
    let abc = ab.cross(ac);
    let acd = ac.cross(ad);
    let adb = ad.cross(ab);

    // Face ABC (opposite to D)
    if abc.dot(ao) > 0.0 {
        // Origin is outside ABC
        simplex.remove(0); // Remove D
        *direction = abc;
        return triangle_case(simplex, direction);
    }

    // Face ACD (opposite to B)
    if acd.dot(ao) > 0.0 {
        // Origin is outside ACD
        simplex.remove(2); // Remove B
        // Reorder: A, C, D -> positions 2, 1, 0
        simplex.points[1] = c;
        simplex.support_a[1] = simplex.support_a[1];
        simplex.support_b[1] = simplex.support_b[1];
        *direction = acd;
        return triangle_case(simplex, direction);
    }

    // Face ADB (opposite to C)
    if adb.dot(ao) > 0.0 {
        // Origin is outside ADB
        simplex.remove(1); // Remove C
        // Reorder to maintain winding
        let temp_point = simplex.points[0];
        let temp_a = simplex.support_a[0];
        let temp_b = simplex.support_b[0];
        simplex.points[0] = simplex.points[1];
        simplex.support_a[0] = simplex.support_a[1];
        simplex.support_b[0] = simplex.support_b[1];
        simplex.points[1] = temp_point;
        simplex.support_a[1] = temp_a;
        simplex.support_b[1] = temp_b;
        *direction = adb;
        return triangle_case(simplex, direction);
    }

    // Origin is inside the tetrahedron
    true
}

/// Computes closest points when shapes are separated
fn compute_closest_points(simplex: &Simplex, _direction: Vec3) -> GjkResult {
    match simplex.size {
        1 => {
            // Single point
            GjkResult::Separated {
                closest_a: simplex.support_a[0],
                closest_b: simplex.support_b[0],
                distance: simplex.points[0].length(),
            }
        }
        2 => {
            // Line segment - find closest point to origin
            let a = simplex.points[1];
            let b = simplex.points[0];
            let ab = b - a;
            let t = (-a).dot(ab) / ab.length_squared();
            let t = t.clamp(0.0, 1.0);

            let closest_a = simplex.support_a[1].lerp(simplex.support_a[0], t);
            let closest_b = simplex.support_b[1].lerp(simplex.support_b[0], t);
            let closest = a.lerp(b, t);

            GjkResult::Separated {
                closest_a,
                closest_b,
                distance: closest.length(),
            }
        }
        3 => {
            // Triangle - find closest point to origin
            let (closest, bary) = closest_point_on_triangle(
                simplex.points[0],
                simplex.points[1],
                simplex.points[2],
            );

            let closest_a = simplex.support_a[0] * bary.x
                + simplex.support_a[1] * bary.y
                + simplex.support_a[2] * bary.z;
            let closest_b = simplex.support_b[0] * bary.x
                + simplex.support_b[1] * bary.y
                + simplex.support_b[2] * bary.z;

            GjkResult::Separated {
                closest_a,
                closest_b,
                distance: closest.length(),
            }
        }
        _ => {
            // Should not happen for separated shapes
            GjkResult::Separated {
                closest_a: Vec3::ZERO,
                closest_b: Vec3::ZERO,
                distance: 0.0,
            }
        }
    }
}

/// Finds the closest point on a triangle to the origin
/// Returns (closest_point, barycentric_coordinates)
fn closest_point_on_triangle(a: Vec3, b: Vec3, c: Vec3) -> (Vec3, Vec3) {
    let ab = b - a;
    let ac = c - a;
    let ao = -a;

    let d1 = ab.dot(ao);
    let d2 = ac.dot(ao);

    // Vertex region A
    if d1 <= 0.0 && d2 <= 0.0 {
        return (a, Vec3::new(1.0, 0.0, 0.0));
    }

    let bo = -b;
    let d3 = ab.dot(bo);
    let d4 = ac.dot(bo);

    // Vertex region B
    if d3 >= 0.0 && d4 <= d3 {
        return (b, Vec3::new(0.0, 1.0, 0.0));
    }

    // Edge region AB
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return (a + ab * v, Vec3::new(1.0 - v, v, 0.0));
    }

    let co = -c;
    let d5 = ab.dot(co);
    let d6 = ac.dot(co);

    // Vertex region C
    if d6 >= 0.0 && d5 <= d6 {
        return (c, Vec3::new(0.0, 0.0, 1.0));
    }

    // Edge region AC
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return (a + ac * w, Vec3::new(1.0 - w, 0.0, w));
    }

    // Edge region BC
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (b + (c - b) * w, Vec3::new(0.0, 1.0 - w, w));
    }

    // Inside triangle
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    (a + ab * v + ac * w, Vec3::new(1.0 - v - w, v, w))
}

/// Simple intersection test using GJK
pub fn intersects(
    shape_a: &Shape,
    transform_a: Transform,
    shape_b: &Shape,
    transform_b: Transform,
) -> bool {
    gjk(shape_a, transform_a, shape_b, transform_b).is_intersecting()
}

/// Computes the distance between two shapes
pub fn distance(
    shape_a: &Shape,
    transform_a: Transform,
    shape_b: &Shape,
    transform_b: Transform,
) -> f32 {
    match gjk(shape_a, transform_a, shape_b, transform_b) {
        GjkResult::Intersecting(_) => 0.0,
        GjkResult::Separated { distance, .. } => distance,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Quat;

    #[test]
    fn test_spheres_intersecting() {
        let sphere = Shape::sphere(1.0);
        let t1 = Transform::from_position(Vec3::ZERO);
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));

        let result = gjk(&sphere, t1, &sphere, t2);
        assert!(result.is_intersecting());
    }

    #[test]
    fn test_spheres_separated() {
        let sphere = Shape::sphere(1.0);
        let t1 = Transform::from_position(Vec3::ZERO);
        let t2 = Transform::from_position(Vec3::new(3.0, 0.0, 0.0));

        let result = gjk(&sphere, t1, &sphere, t2);
        assert!(!result.is_intersecting());

        if let GjkResult::Separated { distance, .. } = result {
            assert!((distance - 1.0).abs() < 0.1); // Distance should be ~1.0
        }
    }

    #[test]
    fn test_box_sphere_intersecting() {
        let box_shape = Shape::cuboid(Vec3::ONE);
        let sphere = Shape::sphere(0.5);

        let t1 = Transform::IDENTITY;
        let t2 = Transform::from_position(Vec3::new(1.2, 0.0, 0.0));

        let result = gjk(&box_shape, t1, &sphere, t2);
        assert!(result.is_intersecting());
    }

    #[test]
    fn test_boxes_separated() {
        let box_shape = Shape::cuboid(Vec3::ONE);

        let t1 = Transform::IDENTITY;
        let t2 = Transform::from_position(Vec3::new(3.0, 0.0, 0.0));

        let result = gjk(&box_shape, t1, &box_shape, t2);
        assert!(!result.is_intersecting());
    }

    #[test]
    fn test_rotated_boxes() {
        let box_shape = Shape::cuboid(Vec3::new(2.0, 0.5, 0.5));

        let t1 = Transform::IDENTITY;
        let t2 = Transform::new(
            Vec3::new(2.0, 0.0, 0.0),
            Quat::from_axis_angle(Vec3::Z, std::f32::consts::FRAC_PI_4),
        );

        // Should intersect when rotated
        let result = gjk(&box_shape, t1, &box_shape, t2);
        assert!(result.is_intersecting());
    }
}
