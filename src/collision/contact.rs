use crate::math::Vec3;

/// Maximum number of contact points in a manifold
pub const MAX_CONTACT_POINTS: usize = 4;

/// A single contact point between two bodies
#[derive(Debug, Clone, Copy)]
pub struct ContactPoint {
    /// Contact point on body A in world space
    pub point_a: Vec3,
    /// Contact point on body B in world space
    pub point_b: Vec3,
    /// Contact point in local space of body A
    pub local_point_a: Vec3,
    /// Contact point in local space of body B
    pub local_point_b: Vec3,
    /// Contact normal (pointing from B to A)
    pub normal: Vec3,
    /// Penetration depth (positive when overlapping)
    pub depth: f32,
    /// Accumulated normal impulse (for warm starting)
    pub normal_impulse: f32,
    /// Accumulated tangent impulse 1 (for warm starting)
    pub tangent_impulse_1: f32,
    /// Accumulated tangent impulse 2 (for warm starting)
    pub tangent_impulse_2: f32,
}

impl ContactPoint {
    /// Creates a new contact point
    pub fn new(
        point_a: Vec3,
        point_b: Vec3,
        local_point_a: Vec3,
        local_point_b: Vec3,
        normal: Vec3,
        depth: f32,
    ) -> Self {
        Self {
            point_a,
            point_b,
            local_point_a,
            local_point_b,
            normal,
            depth,
            normal_impulse: 0.0,
            tangent_impulse_1: 0.0,
            tangent_impulse_2: 0.0,
        }
    }

    /// Returns the midpoint of the contact
    pub fn midpoint(&self) -> Vec3 {
        (self.point_a + self.point_b) * 0.5
    }
}

/// A contact manifold storing contact points between two bodies
#[derive(Debug, Clone)]
pub struct ContactManifold {
    /// Body handle A
    pub body_a: BodyHandle,
    /// Body handle B
    pub body_b: BodyHandle,
    /// Contact points
    pub points: [Option<ContactPoint>; MAX_CONTACT_POINTS],
    /// Number of active contact points
    pub num_points: usize,
    /// Shared contact normal (average of point normals)
    pub normal: Vec3,
    /// Combined friction coefficient
    pub friction: f32,
    /// Combined restitution coefficient
    pub restitution: f32,
}

/// A handle to a body in the physics world
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BodyHandle(pub u32);

impl BodyHandle {
    /// Invalid/null body handle
    pub const INVALID: Self = Self(u32::MAX);

    /// Creates a new body handle
    pub fn new(index: u32) -> Self {
        Self(index)
    }

    /// Returns the index of this handle
    pub fn index(self) -> usize {
        self.0 as usize
    }

    /// Returns true if this handle is valid
    pub fn is_valid(self) -> bool {
        self != Self::INVALID
    }
}

impl Default for BodyHandle {
    fn default() -> Self {
        Self::INVALID
    }
}

impl ContactManifold {
    /// Creates a new empty contact manifold
    pub fn new(body_a: BodyHandle, body_b: BodyHandle) -> Self {
        Self {
            body_a,
            body_b,
            points: [None; MAX_CONTACT_POINTS],
            num_points: 0,
            normal: Vec3::ZERO,
            friction: 0.5,
            restitution: 0.0,
        }
    }

    /// Clears all contact points
    pub fn clear(&mut self) {
        self.points = [None; MAX_CONTACT_POINTS];
        self.num_points = 0;
    }

    /// Adds a contact point to the manifold
    pub fn add_point(&mut self, point: ContactPoint) {
        // Check for duplicate points
        for existing in self.points.iter().flatten() {
            let dist = existing.local_point_a.distance_squared(point.local_point_a);
            if dist < 0.01 {
                return; // Point already exists
            }
        }

        if self.num_points < MAX_CONTACT_POINTS {
            // Find empty slot
            for slot in self.points.iter_mut() {
                if slot.is_none() {
                    *slot = Some(point);
                    self.num_points += 1;
                    self.update_normal();
                    return;
                }
            }
        } else {
            // Manifold is full, replace the point with smallest depth
            let replace_idx = self
                .points
                .iter()
                .enumerate()
                .filter_map(|(i, p)| p.map(|p| (i, p.depth)))
                .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                .map(|(i, _)| i);

            if let Some(idx) = replace_idx {
                if point.depth > self.points[idx].unwrap().depth {
                    self.points[idx] = Some(point);
                    self.update_normal();
                }
            }
        }
    }

    /// Updates the shared normal from point normals
    fn update_normal(&mut self) {
        let sum: Vec3 = self
            .points
            .iter()
            .flatten()
            .map(|p| p.normal)
            .fold(Vec3::ZERO, |a, b| a + b);

        if sum.length_squared() > 0.0001 {
            self.normal = sum.normalize();
        }
    }

    /// Iterates over contact points
    pub fn iter(&self) -> impl Iterator<Item = &ContactPoint> {
        self.points.iter().flatten()
    }

    /// Iterates mutably over contact points
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut ContactPoint> {
        self.points.iter_mut().flatten()
    }

    /// Returns the number of active contact points
    pub fn len(&self) -> usize {
        self.num_points
    }

    /// Returns true if the manifold has no points
    pub fn is_empty(&self) -> bool {
        self.num_points == 0
    }

    /// Refreshes the manifold by removing invalid contacts
    /// Returns the number of points removed
    pub fn refresh(
        &mut self,
        position_a: Vec3,
        rotation_a: crate::math::Quat,
        position_b: Vec3,
        rotation_b: crate::math::Quat,
    ) -> usize {
        let mut removed = 0;

        for slot in self.points.iter_mut() {
            if let Some(point) = slot {
                // Update world space points
                let new_point_a = rotation_a.rotate_vec(point.local_point_a) + position_a;
                let new_point_b = rotation_b.rotate_vec(point.local_point_b) + position_b;

                // Check if contact is still valid
                let diff = new_point_a - new_point_b;
                let normal_dist = diff.dot(point.normal);
                let tangent_dist = (diff - point.normal * normal_dist).length();

                // Remove if too separated or drifted too much
                // normal_dist > 0 means separated along normal, large abs value means moved far
                if normal_dist > 0.01 || normal_dist.abs() > 0.5 || tangent_dist > 0.1 {
                    *slot = None;
                    removed += 1;
                } else {
                    // Update positions
                    point.point_a = new_point_a;
                    point.point_b = new_point_b;
                    point.depth = -normal_dist;
                }
            }
        }

        self.num_points -= removed;
        removed
    }

    /// Warm starts the manifold by matching points from a previous manifold
    pub fn warm_start(&mut self, old: &ContactManifold) {
        for point in self.points.iter_mut().flatten() {
            // Find matching point in old manifold
            for old_point in old.points.iter().flatten() {
                let dist = point.local_point_a.distance_squared(old_point.local_point_a);
                if dist < 0.01 {
                    // Transfer accumulated impulses
                    point.normal_impulse = old_point.normal_impulse;
                    point.tangent_impulse_1 = old_point.tangent_impulse_1;
                    point.tangent_impulse_2 = old_point.tangent_impulse_2;
                    break;
                }
            }
        }
    }
}

/// A collision pair identifier
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CollisionPair {
    /// First body (always has smaller handle)
    pub body_a: BodyHandle,
    /// Second body (always has larger handle)
    pub body_b: BodyHandle,
}

impl CollisionPair {
    /// Creates a new collision pair, ensuring consistent ordering
    pub fn new(a: BodyHandle, b: BodyHandle) -> Self {
        if a.0 <= b.0 {
            Self { body_a: a, body_b: b }
        } else {
            Self { body_a: b, body_b: a }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Quat;

    #[test]
    fn test_manifold_add_points() {
        let mut manifold = ContactManifold::new(BodyHandle::new(0), BodyHandle::new(1));

        let point = ContactPoint::new(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.9, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(-0.1, 0.0, 0.0),
            Vec3::X,
            0.1,
        );

        manifold.add_point(point);
        assert_eq!(manifold.len(), 1);
    }

    #[test]
    fn test_manifold_max_points() {
        let mut manifold = ContactManifold::new(BodyHandle::new(0), BodyHandle::new(1));

        for i in 0..6 {
            let point = ContactPoint::new(
                Vec3::new(i as f32, 0.0, 0.0),
                Vec3::new(i as f32 - 0.1, 0.0, 0.0),
                Vec3::new(i as f32, 0.0, 0.0),
                Vec3::new(-0.1, 0.0, 0.0),
                Vec3::X,
                0.1 * i as f32,
            );
            manifold.add_point(point);
        }

        assert!(manifold.len() <= MAX_CONTACT_POINTS);
    }

    #[test]
    fn test_manifold_refresh() {
        let mut manifold = ContactManifold::new(BodyHandle::new(0), BodyHandle::new(1));

        // Both bodies at origin, contact at (0.5, 0, 0)
        let point = ContactPoint::new(
            Vec3::new(0.5, 0.0, 0.0),   // point_a world
            Vec3::new(0.5, 0.0, 0.0),   // point_b world
            Vec3::new(0.5, 0.0, 0.0),   // local_point_a
            Vec3::new(0.5, 0.0, 0.0),   // local_point_b
            Vec3::X,
            0.1,
        );

        manifold.add_point(point);
        assert_eq!(manifold.len(), 1);

        // Refresh with same transforms - both bodies at origin
        let removed = manifold.refresh(Vec3::ZERO, Quat::IDENTITY, Vec3::ZERO, Quat::IDENTITY);
        assert_eq!(removed, 0, "Should not remove point when bodies haven't moved");
        assert_eq!(manifold.len(), 1);

        // Move body B far away - contact should be invalidated
        let removed = manifold.refresh(Vec3::ZERO, Quat::IDENTITY, Vec3::new(10.0, 0.0, 0.0), Quat::IDENTITY);
        assert_eq!(removed, 1, "Should remove point when body moved far away");
        assert_eq!(manifold.len(), 0);
    }

    #[test]
    fn test_collision_pair_ordering() {
        let pair1 = CollisionPair::new(BodyHandle::new(1), BodyHandle::new(2));
        let pair2 = CollisionPair::new(BodyHandle::new(2), BodyHandle::new(1));

        assert_eq!(pair1, pair2);
        assert_eq!(pair1.body_a.0, 1);
        assert_eq!(pair1.body_b.0, 2);
    }
}
