//! # RustPhy
//!
//! A sophisticated 3D rigid body physics engine written in Rust.
//!
//! ## Features
//!
//! - **Rigid Body Dynamics**: Full 3D rigid body simulation with linear and angular motion
//! - **Collision Detection**: GJK + EPA algorithms for accurate collision detection
//! - **Collision Shapes**: Sphere, Box, and Capsule primitives
//! - **Broad Phase**: Bounding Volume Hierarchy (BVH) for efficient spatial queries
//! - **Constraint Solver**: Projected Gauss-Seidel (PGS) iterative solver
//! - **Contact Caching**: Warm starting for stable simulations
//!
//! ## Quick Start
//!
//! ```rust
//! use rustphy::prelude::*;
//!
//! // Create a physics world
//! let mut world = World::default();
//! world.set_gravity(Vec3::new(0.0, -9.81, 0.0));
//!
//! // Create a static floor
//! let floor = world.create_body(RigidBodyDesc::fixed().with_position(Vec3::ZERO));
//! world.attach_collider(floor, Shape::cuboid(Vec3::new(10.0, 0.5, 10.0)), 1.0);
//!
//! // Create a dynamic ball
//! let ball = world.create_body(
//!     RigidBodyDesc::dynamic()
//!         .with_position(Vec3::new(0.0, 5.0, 0.0))
//!         .with_mass(1.0),
//! );
//! world.attach_collider(ball, Shape::sphere(0.5), 1.0);
//!
//! // Simulation loop
//! let dt = 1.0 / 60.0;
//! for _ in 0..600 {
//!     world.step(dt);
//!     let pos = world.body_position(ball);
//!     println!("Ball position: {:?}", pos);
//! }
//! ```

pub mod collision;
pub mod constraints;
pub mod dynamics;
pub mod geometry;
pub mod math;
pub mod solver;
mod world;

pub use world::{RayCastHit, World, WorldConfig};

/// Prelude module for convenient imports
pub mod prelude {
    pub use crate::collision::{BodyHandle, ContactManifold, ContactPoint};
    pub use crate::dynamics::{BodyType, RigidBody, RigidBodyDesc};
    pub use crate::geometry::{Aabb, BoxShape, Capsule, MassProperties, Shape, ShapeType, Sphere};
    pub use crate::math::{Isometry, Mat3, Quat, Transform, Vec3};
    pub use crate::solver::{PgsSolver, SolverConfig};
    pub use crate::world::{RayCastHit, World, WorldConfig};
}
