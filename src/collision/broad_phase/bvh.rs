use crate::collision::contact::BodyHandle;
use crate::geometry::Aabb;

/// A node in the BVH tree
#[derive(Debug, Clone)]
struct BvhNode {
    /// Bounding box for this node
    aabb: Aabb,
    /// Left child index (or leaf body handle if this is a leaf)
    left: u32,
    /// Right child index (u32::MAX if this is a leaf)
    right: u32,
    /// Parent index
    parent: u32,
    /// Height of the subtree
    height: i32,
}

impl BvhNode {
    fn is_leaf(&self) -> bool {
        self.right == u32::MAX
    }
}

/// A Bounding Volume Hierarchy for broad-phase collision detection
#[derive(Debug)]
pub struct Bvh {
    nodes: Vec<BvhNode>,
    root: Option<u32>,
    /// Maps body handles to their leaf node indices
    body_to_node: Vec<Option<u32>>,
    /// Free node list for reuse
    free_list: Vec<u32>,
    /// AABB margin for fat AABBs
    margin: f32,
}

impl Default for Bvh {
    fn default() -> Self {
        Self::new()
    }
}

impl Bvh {
    /// Creates a new empty BVH
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            root: None,
            body_to_node: Vec::new(),
            free_list: Vec::new(),
            margin: 0.1,
        }
    }

    /// Creates a BVH with a specified margin
    pub fn with_margin(margin: f32) -> Self {
        Self {
            margin,
            ..Self::new()
        }
    }

    /// Inserts a body with its AABB into the BVH
    pub fn insert(&mut self, handle: BodyHandle, aabb: Aabb) {
        let fat_aabb = aabb.expand(self.margin);
        let leaf = self.allocate_leaf(handle, fat_aabb);

        // Ensure body_to_node is large enough
        let index = handle.index();
        if index >= self.body_to_node.len() {
            self.body_to_node.resize(index + 1, None);
        }
        self.body_to_node[index] = Some(leaf);

        if self.root.is_none() {
            self.root = Some(leaf);
            return;
        }

        // Find best sibling
        let sibling = self.find_best_sibling(leaf);

        // Create new parent
        let old_parent = self.nodes[sibling as usize].parent;
        let new_parent = self.allocate_internal(
            self.nodes[leaf as usize].aabb.union(self.nodes[sibling as usize].aabb),
        );

        self.nodes[new_parent as usize].parent = old_parent;
        self.nodes[new_parent as usize].left = sibling;
        self.nodes[new_parent as usize].right = leaf;
        self.nodes[sibling as usize].parent = new_parent;
        self.nodes[leaf as usize].parent = new_parent;

        if old_parent == u32::MAX {
            self.root = Some(new_parent);
        } else {
            let old_parent_node = &mut self.nodes[old_parent as usize];
            if old_parent_node.left == sibling {
                old_parent_node.left = new_parent;
            } else {
                old_parent_node.right = new_parent;
            }
        }

        // Refit ancestors
        self.refit(new_parent);
    }

    /// Removes a body from the BVH
    pub fn remove(&mut self, handle: BodyHandle) {
        let index = handle.index();
        if index >= self.body_to_node.len() {
            return;
        }

        let leaf = match self.body_to_node[index] {
            Some(l) => l,
            None => return,
        };
        self.body_to_node[index] = None;

        if Some(leaf) == self.root {
            self.root = None;
            self.free_node(leaf);
            return;
        }

        let parent = self.nodes[leaf as usize].parent;
        let grandparent = self.nodes[parent as usize].parent;
        let sibling = if self.nodes[parent as usize].left == leaf {
            self.nodes[parent as usize].right
        } else {
            self.nodes[parent as usize].left
        };

        if grandparent == u32::MAX {
            self.root = Some(sibling);
            self.nodes[sibling as usize].parent = u32::MAX;
        } else {
            let grandparent_node = &mut self.nodes[grandparent as usize];
            if grandparent_node.left == parent {
                grandparent_node.left = sibling;
            } else {
                grandparent_node.right = sibling;
            }
            self.nodes[sibling as usize].parent = grandparent;
            self.refit(grandparent);
        }

        self.free_node(leaf);
        self.free_node(parent);
    }

    /// Updates a body's AABB
    pub fn update(&mut self, handle: BodyHandle, aabb: Aabb) -> bool {
        let index = handle.index();
        if index >= self.body_to_node.len() {
            return false;
        }

        let leaf = match self.body_to_node[index] {
            Some(l) => l,
            None => return false,
        };

        let fat_aabb = &self.nodes[leaf as usize].aabb;

        // Check if the fat AABB still contains the new AABB
        if fat_aabb.contains_aabb(aabb) {
            return false; // No update needed
        }

        // Remove and reinsert
        self.remove(handle);
        self.insert(handle, aabb);
        true
    }

    /// Queries all pairs of potentially colliding bodies
    pub fn query_pairs(&self) -> Vec<(BodyHandle, BodyHandle)> {
        let mut pairs = Vec::new();

        if self.root.is_none() {
            return pairs;
        }

        // Collect all leaf nodes
        let mut leaves = Vec::new();
        self.collect_leaves(self.root.unwrap(), &mut leaves);

        // Test each leaf against the tree
        for &leaf in &leaves {
            let handle = BodyHandle::new(self.nodes[leaf as usize].left);
            let aabb = self.nodes[leaf as usize].aabb;

            self.query_aabb_internal(self.root.unwrap(), aabb, |other| {
                if other != handle && handle.0 < other.0 {
                    pairs.push((handle, other));
                }
            });
        }

        pairs
    }

    /// Queries bodies whose AABBs intersect the given AABB
    pub fn query_aabb(&self, aabb: Aabb, callback: impl FnMut(BodyHandle)) {
        if let Some(root) = self.root {
            self.query_aabb_internal(root, aabb, callback);
        }
    }

    fn query_aabb_internal(&self, node: u32, aabb: Aabb, mut callback: impl FnMut(BodyHandle)) {
        let mut stack = vec![node];

        while let Some(current) = stack.pop() {
            let n = &self.nodes[current as usize];

            if !n.aabb.intersects(aabb) {
                continue;
            }

            if n.is_leaf() {
                callback(BodyHandle::new(n.left));
            } else {
                stack.push(n.left);
                stack.push(n.right);
            }
        }
    }

    /// Queries bodies that a ray might intersect
    pub fn query_ray(
        &self,
        origin: crate::math::Vec3,
        direction: crate::math::Vec3,
        max_distance: f32,
    ) -> Vec<BodyHandle> {
        let mut results = Vec::new();

        if self.root.is_none() {
            return results;
        }

        let mut stack = vec![self.root.unwrap()];

        while let Some(current) = stack.pop() {
            let n = &self.nodes[current as usize];

            if let Some((t_min, _)) = n.aabb.ray_intersection(origin, direction) {
                if t_min > max_distance {
                    continue;
                }

                if n.is_leaf() {
                    results.push(BodyHandle::new(n.left));
                } else {
                    stack.push(n.left);
                    stack.push(n.right);
                }
            }
        }

        results
    }

    fn collect_leaves(&self, node: u32, leaves: &mut Vec<u32>) {
        let n = &self.nodes[node as usize];
        if n.is_leaf() {
            leaves.push(node);
        } else {
            self.collect_leaves(n.left, leaves);
            self.collect_leaves(n.right, leaves);
        }
    }

    fn allocate_leaf(&mut self, handle: BodyHandle, aabb: Aabb) -> u32 {
        let node = BvhNode {
            aabb,
            left: handle.0,
            right: u32::MAX, // Marks as leaf
            parent: u32::MAX,
            height: 0,
        };

        self.allocate_node(node)
    }

    fn allocate_internal(&mut self, aabb: Aabb) -> u32 {
        let node = BvhNode {
            aabb,
            left: u32::MAX,
            right: u32::MAX,
            parent: u32::MAX,
            height: 0,
        };

        self.allocate_node(node)
    }

    fn allocate_node(&mut self, node: BvhNode) -> u32 {
        if let Some(index) = self.free_list.pop() {
            self.nodes[index as usize] = node;
            index
        } else {
            let index = self.nodes.len() as u32;
            self.nodes.push(node);
            index
        }
    }

    fn free_node(&mut self, index: u32) {
        self.free_list.push(index);
    }

    fn find_best_sibling(&self, leaf: u32) -> u32 {
        let leaf_aabb = self.nodes[leaf as usize].aabb;
        let mut best = self.root.unwrap();
        let mut best_cost = leaf_aabb.union(self.nodes[best as usize].aabb).surface_area();

        let mut stack = vec![self.root.unwrap()];

        while let Some(current) = stack.pop() {
            let n = &self.nodes[current as usize];

            let combined = leaf_aabb.union(n.aabb);
            let combined_cost = combined.surface_area();

            if combined_cost < best_cost {
                best = current;
                best_cost = combined_cost;
            }

            if !n.is_leaf() {
                // Estimate cost of descending
                let inherited_cost = combined_cost - n.aabb.surface_area();
                let left_cost = leaf_aabb.union(self.nodes[n.left as usize].aabb).surface_area()
                    + inherited_cost;
                let right_cost = leaf_aabb.union(self.nodes[n.right as usize].aabb).surface_area()
                    + inherited_cost;

                if left_cost < best_cost || right_cost < best_cost {
                    stack.push(n.left);
                    stack.push(n.right);
                }
            }
        }

        best
    }

    fn refit(&mut self, start: u32) {
        let mut current = start;

        while current != u32::MAX {
            let n = &self.nodes[current as usize];
            if n.is_leaf() {
                current = n.parent;
                continue;
            }

            let left = n.left;
            let right = n.right;
            let left_aabb = self.nodes[left as usize].aabb;
            let right_aabb = self.nodes[right as usize].aabb;
            let left_height = self.nodes[left as usize].height;
            let right_height = self.nodes[right as usize].height;

            self.nodes[current as usize].aabb = left_aabb.union(right_aabb);
            self.nodes[current as usize].height = 1 + left_height.max(right_height);

            current = self.nodes[current as usize].parent;
        }
    }

    /// Returns the number of bodies in the BVH
    pub fn len(&self) -> usize {
        self.body_to_node.iter().filter(|x| x.is_some()).count()
    }

    /// Returns true if the BVH is empty
    pub fn is_empty(&self) -> bool {
        self.root.is_none()
    }

    /// Clears all bodies from the BVH
    pub fn clear(&mut self) {
        self.nodes.clear();
        self.root = None;
        self.body_to_node.clear();
        self.free_list.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Vec3;

    #[test]
    fn test_insert_and_query() {
        let mut bvh = Bvh::new();

        let aabb1 = Aabb::new(Vec3::ZERO, Vec3::ONE);
        let aabb2 = Aabb::new(Vec3::new(0.5, 0.0, 0.0), Vec3::new(1.5, 1.0, 1.0));
        let aabb3 = Aabb::new(Vec3::new(5.0, 0.0, 0.0), Vec3::new(6.0, 1.0, 1.0));

        bvh.insert(BodyHandle::new(0), aabb1);
        bvh.insert(BodyHandle::new(1), aabb2);
        bvh.insert(BodyHandle::new(2), aabb3);

        let pairs = bvh.query_pairs();

        // Bodies 0 and 1 should overlap
        assert!(pairs.contains(&(BodyHandle::new(0), BodyHandle::new(1))));

        // Body 2 shouldn't overlap with others
        assert!(!pairs.contains(&(BodyHandle::new(0), BodyHandle::new(2))));
        assert!(!pairs.contains(&(BodyHandle::new(1), BodyHandle::new(2))));
    }

    #[test]
    fn test_remove() {
        let mut bvh = Bvh::new();

        let aabb1 = Aabb::new(Vec3::ZERO, Vec3::ONE);
        let aabb2 = Aabb::new(Vec3::new(0.5, 0.0, 0.0), Vec3::new(1.5, 1.0, 1.0));

        bvh.insert(BodyHandle::new(0), aabb1);
        bvh.insert(BodyHandle::new(1), aabb2);

        assert_eq!(bvh.len(), 2);

        bvh.remove(BodyHandle::new(0));
        assert_eq!(bvh.len(), 1);

        let pairs = bvh.query_pairs();
        assert!(pairs.is_empty());
    }

    #[test]
    fn test_update() {
        let mut bvh = Bvh::new();

        let aabb1 = Aabb::new(Vec3::ZERO, Vec3::ONE);
        let aabb2 = Aabb::new(Vec3::new(5.0, 0.0, 0.0), Vec3::new(6.0, 1.0, 1.0));

        bvh.insert(BodyHandle::new(0), aabb1);
        bvh.insert(BodyHandle::new(1), aabb2);

        // Initially no overlap
        let pairs = bvh.query_pairs();
        assert!(!pairs.contains(&(BodyHandle::new(0), BodyHandle::new(1))));

        // Move body 1 to overlap with body 0
        let new_aabb2 = Aabb::new(Vec3::new(0.5, 0.0, 0.0), Vec3::new(1.5, 1.0, 1.0));
        bvh.update(BodyHandle::new(1), new_aabb2);

        let pairs = bvh.query_pairs();
        assert!(pairs.contains(&(BodyHandle::new(0), BodyHandle::new(1))));
    }

    #[test]
    fn test_ray_query() {
        let mut bvh = Bvh::new();

        bvh.insert(BodyHandle::new(0), Aabb::new(Vec3::ZERO, Vec3::ONE));
        bvh.insert(
            BodyHandle::new(1),
            Aabb::new(Vec3::new(5.0, 0.0, 0.0), Vec3::new(6.0, 1.0, 1.0)),
        );

        let hits = bvh.query_ray(Vec3::new(-1.0, 0.5, 0.5), Vec3::X, 100.0);

        assert!(hits.contains(&BodyHandle::new(0)));
        assert!(hits.contains(&BodyHandle::new(1)));
    }
}
