use crate::collision::{ContactManifold, ContactPoint};
use crate::constraints::{warm_start, ContactVelocityConstraint};
use crate::dynamics::RigidBody;

/// Configuration for the constraint solver
#[derive(Debug, Clone, Copy)]
pub struct SolverConfig {
    /// Number of velocity solver iterations
    pub velocity_iterations: usize,
    /// Number of position solver iterations
    pub position_iterations: usize,
    /// Warm starting coefficient (0-1)
    pub warm_start_coefficient: f32,
    /// Baumgarte stabilization coefficient
    pub baumgarte: f32,
    /// Allowed penetration slop
    pub slop: f32,
}

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            velocity_iterations: 8,    // Iterations for convergence
            position_iterations: 4,    // More position correction iterations
            warm_start_coefficient: 0.8, // Warm starting
            baumgarte: 0.5,            // Strong position correction (position solver only)
            slop: 0.005,               // Small slop
        }
    }
}

/// Projected Gauss-Seidel constraint solver
pub struct PgsSolver {
    config: SolverConfig,
    velocity_constraints: Vec<(usize, usize, Vec<ContactVelocityConstraint>)>,
}

impl Default for PgsSolver {
    fn default() -> Self {
        Self::new(SolverConfig::default())
    }
}

impl PgsSolver {
    /// Creates a new PGS solver with the given configuration
    pub fn new(config: SolverConfig) -> Self {
        Self {
            config,
            velocity_constraints: Vec::new(),
        }
    }

    /// Prepares velocity constraints from contact manifolds
    pub fn prepare(
        &mut self,
        manifolds: &[ContactManifold],
        bodies: &[RigidBody],
        dt: f32,
    ) {
        self.velocity_constraints.clear();

        for manifold in manifolds {
            let body_a_idx = manifold.body_a.index();
            let body_b_idx = manifold.body_b.index();

            if body_a_idx >= bodies.len() || body_b_idx >= bodies.len() {
                continue;
            }

            let body_a = &bodies[body_a_idx];
            let body_b = &bodies[body_b_idx];

            // Skip if both bodies are static
            if !body_a.is_dynamic() && !body_b.is_dynamic() {
                continue;
            }

            let restitution = (body_a.restitution + body_b.restitution) * 0.5;

            let mut constraints = Vec::new();

            for contact in manifold.iter() {
                let constraint =
                    ContactVelocityConstraint::new(contact, body_a, body_b, restitution, dt);
                constraints.push(constraint);
            }

            if !constraints.is_empty() {
                self.velocity_constraints.push((body_a_idx, body_b_idx, constraints));
            }
        }
    }

    /// Warm starts the solver using previously accumulated impulses
    pub fn warm_start(&self, bodies: &mut [RigidBody]) {
        for (body_a_idx, body_b_idx, constraints) in &self.velocity_constraints {
            for constraint in constraints {
                let scaled_constraint = ContactVelocityConstraint {
                    normal_impulse: constraint.normal_impulse * self.config.warm_start_coefficient,
                    tangent1_impulse: constraint.tangent1_impulse * self.config.warm_start_coefficient,
                    tangent2_impulse: constraint.tangent2_impulse * self.config.warm_start_coefficient,
                    ..*constraint
                };

                let (bodies_before, bodies_after) = bodies.split_at_mut(*body_a_idx.max(body_b_idx));
                let (body_a, body_b) = if body_a_idx < body_b_idx {
                    (&mut bodies_before[*body_a_idx], &mut bodies_after[0])
                } else {
                    (&mut bodies_after[0], &mut bodies_before[*body_b_idx])
                };

                warm_start(&scaled_constraint, body_a, body_b);
            }
        }
    }

    /// Solves velocity constraints
    pub fn solve_velocity(&mut self, bodies: &mut [RigidBody]) {
        for _ in 0..self.config.velocity_iterations {
            for (body_a_idx, body_b_idx, constraints) in &mut self.velocity_constraints {
                let friction = {
                    let body_a = &bodies[*body_a_idx];
                    let body_b = &bodies[*body_b_idx];
                    (body_a.friction * body_b.friction).sqrt()
                };

                for constraint in constraints.iter_mut() {
                    // Get mutable references to both bodies
                    let (body_a, body_b) = get_two_mut(bodies, *body_a_idx, *body_b_idx);

                    // Solve normal constraint
                    constraint.solve_normal(body_a, body_b);

                    // Solve friction constraints
                    constraint.solve_friction(body_a, body_b, friction);
                }
            }
        }
    }

    /// Stores accumulated impulses back to contact manifolds
    pub fn store_impulses(&self, manifolds: &mut [ContactManifold]) {
        for (i, (_, _, constraints)) in self.velocity_constraints.iter().enumerate() {
            if i < manifolds.len() {
                let manifold = &mut manifolds[i];
                for (j, constraint) in constraints.iter().enumerate() {
                    if let Some(contact) = manifold.points.get_mut(j).and_then(|p| p.as_mut()) {
                        constraint.store_impulses(contact);
                    }
                }
            }
        }
    }

    /// Returns the solver configuration
    pub fn config(&self) -> &SolverConfig {
        &self.config
    }

    /// Sets the solver configuration
    pub fn set_config(&mut self, config: SolverConfig) {
        self.config = config;
    }
}

/// Gets mutable references to two elements at different indices
fn get_two_mut(slice: &mut [RigidBody], a: usize, b: usize) -> (&mut RigidBody, &mut RigidBody) {
    assert!(a != b);
    if a < b {
        let (left, right) = slice.split_at_mut(b);
        (&mut left[a], &mut right[0])
    } else {
        let (left, right) = slice.split_at_mut(a);
        (&mut right[0], &mut left[b])
    }
}

/// Solves position constraints to resolve penetration
pub fn solve_position_constraints(
    manifolds: &[ContactManifold],
    bodies: &mut [RigidBody],
    config: &SolverConfig,
) -> bool {
    let mut solved = true;

    for _ in 0..config.position_iterations {
        let mut max_penetration = 0.0f32;

        for manifold in manifolds {
            let body_a_idx = manifold.body_a.index();
            let body_b_idx = manifold.body_b.index();

            if body_a_idx >= bodies.len() || body_b_idx >= bodies.len() {
                continue;
            }

            // Skip if both bodies are static
            if !bodies[body_a_idx].is_dynamic() && !bodies[body_b_idx].is_dynamic() {
                continue;
            }

            for contact in manifold.iter() {
                let penetration = contact.depth - config.slop;
                max_penetration = max_penetration.max(penetration);

                if penetration > 0.0 {
                    let (body_a, body_b) = get_two_mut(bodies, body_a_idx, body_b_idx);
                    solve_single_position_constraint(contact, body_a, body_b, config);
                }
            }
        }

        if max_penetration < config.slop {
            break;
        }

        if max_penetration > config.slop * 2.0 {
            solved = false;
        }
    }

    solved
}

/// Solves a single position constraint
fn solve_single_position_constraint(
    contact: &ContactPoint,
    body_a: &mut RigidBody,
    body_b: &mut RigidBody,
    config: &SolverConfig,
) {
    let penetration = (contact.depth - config.slop).max(0.0);
    if penetration <= 0.0 {
        return;
    }

    let point = (contact.point_a + contact.point_b) * 0.5;
    let normal = contact.normal;

    let r_a = point - body_a.position;
    let r_b = point - body_b.position;

    // Compute effective mass
    let rn_a = r_a.cross(normal);
    let rn_b = r_b.cross(normal);

    let k = body_a.inv_mass
        + body_b.inv_mass
        + rn_a.dot(body_a.inv_inertia_world * rn_a)
        + rn_b.dot(body_b.inv_inertia_world * rn_b);

    if k <= 0.0 {
        return;
    }

    let correction = penetration * config.baumgarte / k;

    // Apply position correction
    // Normal points from B to A, so:
    // - B should move in +normal direction (away from A)
    // - A should move in -normal direction (away from B)
    body_a.position -= normal * (correction * body_a.inv_mass);
    body_b.position += normal * (correction * body_b.inv_mass);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::collision::BodyHandle;
    use crate::dynamics::BodyType;
    use crate::math::Vec3;

    #[test]
    fn test_solver_creation() {
        let solver = PgsSolver::default();
        assert_eq!(solver.config().velocity_iterations, 8);
    }

    #[test]
    fn test_prepare_constraints() {
        let mut solver = PgsSolver::default();

        let mut manifold = ContactManifold::new(BodyHandle::new(0), BodyHandle::new(1));
        let contact = ContactPoint::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, -0.1, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, -0.1, 0.0),
            Vec3::Y,
            0.1,
        );
        manifold.add_point(contact);

        let mut body_a = RigidBody::new().with_mass(1.0);
        let body_b = RigidBody::new().with_type(BodyType::Static);
        body_a.update_world_inertia();

        let bodies = vec![body_a, body_b];
        let manifolds = vec![manifold];

        solver.prepare(&manifolds, &bodies, 1.0 / 60.0);

        assert_eq!(solver.velocity_constraints.len(), 1);
    }
}
