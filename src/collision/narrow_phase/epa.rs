use crate::geometry::Shape;
use crate::math::{Transform, Vec3};

use super::gjk::Simplex;

/// Maximum iterations for EPA algorithm
const EPA_MAX_ITERATIONS: usize = 64;

/// Tolerance for EPA convergence
const EPA_TOLERANCE: f32 = 1e-4;

/// Maximum number of faces in the polytope
const EPA_MAX_FACES: usize = 128;

/// Result of EPA algorithm
#[derive(Debug, Clone, Copy)]
pub struct EpaResult {
    /// Penetration normal (pointing from B to A)
    pub normal: Vec3,
    /// Penetration depth
    pub depth: f32,
    /// Contact point on shape A (in world space)
    pub point_a: Vec3,
    /// Contact point on shape B (in world space)
    pub point_b: Vec3,
}

/// A face of the polytope
#[derive(Debug, Clone, Copy)]
struct Face {
    /// Indices of the three vertices
    indices: [usize; 3],
    /// Face normal (pointing outward)
    normal: Vec3,
    /// Distance from origin to the face
    distance: f32,
}

/// A vertex in the polytope
#[derive(Debug, Clone, Copy)]
struct Vertex {
    /// Point in Minkowski difference space
    point: Vec3,
    /// Support point on shape A
    support_a: Vec3,
    /// Support point on shape B
    support_b: Vec3,
}

/// Performs the EPA algorithm to find penetration depth and contact points
pub fn epa(
    simplex: &Simplex,
    shape_a: &Shape,
    transform_a: Transform,
    shape_b: &Shape,
    transform_b: Transform,
) -> Option<EpaResult> {
    if simplex.size < 4 {
        return None;
    }

    // Initialize polytope from the GJK simplex
    let mut vertices: Vec<Vertex> = Vec::with_capacity(EPA_MAX_FACES);
    let mut faces: Vec<Face> = Vec::with_capacity(EPA_MAX_FACES);

    // Add simplex vertices
    for i in 0..4 {
        vertices.push(Vertex {
            point: simplex.points[i],
            support_a: simplex.support_a[i],
            support_b: simplex.support_b[i],
        });
    }

    // Create initial tetrahedron faces (ensure outward normals)
    let faces_indices = [
        [0, 1, 2],
        [0, 3, 1],
        [0, 2, 3],
        [1, 3, 2],
    ];

    for indices in faces_indices.iter() {
        if let Some(face) = create_face(&vertices, *indices) {
            faces.push(face);
        }
    }

    // Fix winding order to ensure normals point outward
    fix_winding(&vertices, &mut faces);

    for _ in 0..EPA_MAX_ITERATIONS {
        if faces.is_empty() {
            return None;
        }

        // Find the closest face to the origin
        let (closest_idx, closest_face) = faces
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.distance.partial_cmp(&b.distance).unwrap())
            .map(|(i, f)| (i, *f))?;

        // Get support point in the direction of the closest face normal
        let (new_point, sup_a, sup_b) =
            support(shape_a, transform_a, shape_b, transform_b, closest_face.normal);

        // Check for convergence
        let distance = new_point.dot(closest_face.normal);
        if distance - closest_face.distance < EPA_TOLERANCE {
            // Converged - compute contact information
            return Some(compute_contact(&vertices, &closest_face));
        }

        // Add new vertex
        let new_vertex_idx = vertices.len();
        vertices.push(Vertex {
            point: new_point,
            support_a: sup_a,
            support_b: sup_b,
        });

        // Find and remove faces visible from the new point
        let mut edges_to_patch: Vec<(usize, usize)> = Vec::new();
        let mut faces_to_remove: Vec<usize> = Vec::new();

        for (i, face) in faces.iter().enumerate() {
            let v = vertices[face.indices[0]].point;
            if face.normal.dot(new_point - v) > 0.0 {
                // Face is visible from new point, mark for removal
                faces_to_remove.push(i);

                // Collect edges
                for j in 0..3 {
                    let edge = (face.indices[j], face.indices[(j + 1) % 3]);
                    add_or_remove_edge(&mut edges_to_patch, edge);
                }
            }
        }

        // Remove faces (in reverse order to maintain indices)
        faces_to_remove.sort_by(|a, b| b.cmp(a));
        for idx in faces_to_remove {
            faces.swap_remove(idx);
        }

        // Create new faces from the edges
        for (a, b) in edges_to_patch {
            if let Some(face) = create_face(&vertices, [a, b, new_vertex_idx]) {
                faces.push(face);
            }
        }

        if faces.len() > EPA_MAX_FACES {
            // Too many faces, return best result so far
            let closest_face = faces
                .iter()
                .min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap())?;
            return Some(compute_contact(&vertices, closest_face));
        }
    }

    // Max iterations reached
    let closest_face = faces
        .iter()
        .min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap())?;
    Some(compute_contact(&vertices, closest_face))
}

/// Computes support point in Minkowski difference
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

/// Creates a face from three vertex indices
fn create_face(vertices: &[Vertex], indices: [usize; 3]) -> Option<Face> {
    let a = vertices[indices[0]].point;
    let b = vertices[indices[1]].point;
    let c = vertices[indices[2]].point;

    let ab = b - a;
    let ac = c - a;
    let normal = ab.cross(ac);

    let len = normal.length();
    if len < 1e-10 {
        return None; // Degenerate face
    }

    let normal = normal / len;
    let distance = normal.dot(a);

    Some(Face {
        indices,
        normal,
        distance,
    })
}

/// Fixes face winding to ensure normals point outward (away from origin)
fn fix_winding(vertices: &[Vertex], faces: &mut [Face]) {
    // Compute centroid of the polytope
    let centroid: Vec3 = vertices.iter().map(|v| v.point).fold(Vec3::ZERO, |a, b| a + b)
        / vertices.len() as f32;

    for face in faces.iter_mut() {
        let v = vertices[face.indices[0]].point;
        let to_face = v - centroid;

        // If normal points toward centroid, flip it
        if face.normal.dot(to_face) < 0.0 {
            face.normal = -face.normal;
            face.distance = -face.distance;
            // Swap two vertices to fix winding
            face.indices.swap(0, 1);
        }
    }
}

/// Adds an edge to the list, or removes it if it already exists (shared edge)
fn add_or_remove_edge(edges: &mut Vec<(usize, usize)>, edge: (usize, usize)) {
    let reverse = (edge.1, edge.0);

    if let Some(pos) = edges.iter().position(|e| *e == reverse) {
        edges.remove(pos);
    } else {
        edges.push(edge);
    }
}

/// Computes contact information from the closest face
fn compute_contact(vertices: &[Vertex], face: &Face) -> EpaResult {
    let a = vertices[face.indices[0]].point;
    let b = vertices[face.indices[1]].point;
    let c = vertices[face.indices[2]].point;

    // Project origin onto the face to get barycentric coordinates
    let bary = barycentric_coordinates(Vec3::ZERO, a, b, c);

    // Compute contact points using barycentric interpolation
    let point_a = vertices[face.indices[0]].support_a * bary.x
        + vertices[face.indices[1]].support_a * bary.y
        + vertices[face.indices[2]].support_a * bary.z;

    let point_b = vertices[face.indices[0]].support_b * bary.x
        + vertices[face.indices[1]].support_b * bary.y
        + vertices[face.indices[2]].support_b * bary.z;

    EpaResult {
        normal: face.normal,
        depth: face.distance.abs(),
        point_a,
        point_b,
    }
}

/// Computes barycentric coordinates for a point with respect to a triangle
fn barycentric_coordinates(p: Vec3, a: Vec3, b: Vec3, c: Vec3) -> Vec3 {
    let v0 = b - a;
    let v1 = c - a;
    let v2 = p - a;

    let d00 = v0.dot(v0);
    let d01 = v0.dot(v1);
    let d11 = v1.dot(v1);
    let d20 = v2.dot(v0);
    let d21 = v2.dot(v1);

    let denom = d00 * d11 - d01 * d01;

    if denom.abs() < 1e-10 {
        // Degenerate triangle, return equal weights
        return Vec3::new(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    }

    let v = (d11 * d20 - d01 * d21) / denom;
    let w = (d00 * d21 - d01 * d20) / denom;
    let u = 1.0 - v - w;

    Vec3::new(u, v, w)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::collision::narrow_phase::gjk::{gjk, GjkResult};

    #[test]
    fn test_sphere_sphere_penetration() {
        let sphere = Shape::sphere(1.0);
        let t1 = Transform::from_position(Vec3::ZERO);
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));

        let result = gjk(&sphere, t1, &sphere, t2);

        if let GjkResult::Intersecting(simplex) = result {
            let epa_result = epa(&simplex, &sphere, t1, &sphere, t2);
            assert!(epa_result.is_some());

            let contact = epa_result.unwrap();
            // Penetration depth should be approximately 0.5 (2*1.0 - 1.5)
            assert!((contact.depth - 0.5).abs() < 0.1);
            // Normal should point along X axis
            assert!(contact.normal.x.abs() > 0.9);
        } else {
            panic!("Expected intersection");
        }
    }

    #[test]
    fn test_box_box_penetration() {
        let box_shape = Shape::cuboid(Vec3::ONE);
        let t1 = Transform::IDENTITY;
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));

        let result = gjk(&box_shape, t1, &box_shape, t2);

        if let GjkResult::Intersecting(simplex) = result {
            let epa_result = epa(&simplex, &box_shape, t1, &box_shape, t2);
            assert!(epa_result.is_some());

            let contact = epa_result.unwrap();
            // Penetration depth should be approximately 0.5
            assert!((contact.depth - 0.5).abs() < 0.1);
        } else {
            panic!("Expected intersection");
        }
    }
}
