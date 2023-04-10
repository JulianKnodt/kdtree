#![allow(unused)]
type F = f32;

fn dist_sq<const N: usize>(a: &[F; N], b: &[F; N]) -> F {
    let mut d = 0.;
    for i in 0..N {
        let v = a[i] - b[i];
        d += v * v;
    }
    d
}

fn dist<const N: usize>(a: &[F; N], b: &[F; N]) -> F {
    dist_sq(a, b).sqrt()
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct AABB<const N: usize> {
    min: [F; N],
    max: [F; N],
}

impl<const N: usize> AABB<N> {
    const EMPTY: Self = AABB {
        min: [F::INFINITY; N],
        max: [F::NEG_INFINITY; N],
    };
    pub fn add_point(&mut self, p: &[F; N]) {
        for i in 0..N {
            self.min[i] = p[i].min(self.min[i]);
            self.max[i] = p[i].max(self.max[i]);
        }
    }
    pub fn center(&self) -> [F; N] {
        std::array::from_fn(|i| (self.min[i] + self.max[i]) / 2.)
    }
    pub fn extent_length(&self) -> F {
        dist(&self.max, &self.min)
    }
    pub fn to_sphere(&self) -> Sphere<N> {
        Sphere {
            center: self.center(),
            radius: self.extent_length() / 2.,
        }
    }
    pub fn largest_dimension(&self) -> usize {
        (0..N)
            .map(|i| (i, self.max[i] - self.min[i]))
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .unwrap()
            .0
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Sphere<const N: usize> {
    center: [F; N],
    radius: F,
}

impl<const N: usize> Sphere<N> {
    const EMPTY: Self = Self {
        center: [0.; N],
        radius: F::INFINITY,
    };
    // If it overlaps, returns the distance to the sphere
    fn overlaps(&self, pt: &[F; N], rad: F) -> Option<F> {
        let d = dist(&self.center, pt);
        if d <= rad + self.radius {
            Some(d - self.radius)
        } else {
            None
        }
    }
    fn volume(&self) -> F {
        const PI: F = std::f64::consts::PI as F;
        4. / 3. * PI * self.radius * self.radius * self.radius
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
struct KDNode<const N: usize> {
    bounds: Sphere<N>,

    left_child_or_first_point: usize,
    num_points: usize,
}

impl<const N: usize> KDNode<N> {
    const EMPTY: Self = KDNode {
        bounds: Sphere::EMPTY,
        left_child_or_first_point: 0,
        num_points: 0,
    };

    fn right_child(&self) -> usize {
        assert!(!self.is_leaf());
        self.left_child() + 1
    }
    fn is_leaf(&self) -> bool {
        self.num_points > 0
    }
    fn set_left_child(&mut self, left_child: usize) -> (usize, usize) {
        assert!(self.is_leaf());
        let old_first_pt = std::mem::replace(&mut self.left_child_or_first_point, left_child);
        let num_prims = std::mem::take(&mut self.num_points);
        (old_first_pt, num_prims)
    }
    fn left_child(&self) -> usize {
        assert!(!self.is_leaf());
        self.left_child_or_first_point
    }
    fn first_point(&self) -> usize {
        assert!(self.is_leaf());
        self.left_child_or_first_point
    }
    fn set_points(&mut self, first_pt: usize, num_pts: usize) {
        self.left_child_or_first_point = first_pt;
        self.num_points = num_pts;
    }
}

#[derive(Debug, Clone, Default)]
pub struct KDTree<T, Q, const N: usize> {
    nodes: Vec<KDNode<N>>,
    root_node_idx: usize,
    nodes_used: usize,

    points: Vec<[T; N]>,
    data: Vec<Q>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SplitKind {
    /// Split along the middle of each sphere
    Midpoint,
    // TODO add SAH and other volume heuristic
}

impl From<()> for SplitKind {
    fn from((): ()) -> SplitKind {
        SplitKind::Midpoint
    }
}

impl<const N: usize> KDTree<F, (), N> {
    pub fn new(pts: impl Iterator<Item = [F; N]>, split: impl Into<SplitKind>) -> Self {
        let points = pts.collect::<Vec<_>>();
        let size = 2 * points.len() + 1;
        let nodes = vec![KDNode::EMPTY; size];
        let mut s = Self {
            data: vec![(); points.len()],
            nodes,
            root_node_idx: 0,
            nodes_used: 1,
            points,
        };
        s.nodes[0].num_points = s.points.len();
        s.update_node_bounds(0);
        s.subdivide(0, split.into());
        s
    }
    fn update_node_bounds(&mut self, idx: usize) {
        let node = &mut self.nodes[idx];
        let mut aabb = AABB::EMPTY;
        let fp = node.first_point();
        for p in &self.points[fp..fp + node.num_points] {
            aabb.add_point(p);
        }
        let mut sphere = aabb.to_sphere();
        sphere.radius += 1e-7;
        node.bounds = sphere;
    }
    fn midpoint_split(&self, node: &KDNode<N>) -> (usize, F) {
        let mut aabb = AABB::EMPTY;
        let fp = node.first_point();
        for p in &self.points[fp..fp + node.num_points] {
            aabb.add_point(p);
        }
        let axis = aabb.largest_dimension();
        (axis, aabb.center()[axis])
    }
    fn subdivide(&mut self, idx: usize, split_kind: SplitKind) {
        let node = &self.nodes[idx];
        // TODO here can use a different amount of points so it will be faster?
        if node.num_points <= 16 {
            return;
        }
        let (axis, split_val) = match split_kind {
            SplitKind::Midpoint => self.midpoint_split(node),
        };
        /*
        let (cost,axis,split_pos) = self.sah_linplace_split_binned::<2048>(node);
        let curr_cost = node.aabb().volume() * (node.num_prims as F)
        if cost >= curr_cost {
          return;
        }
        */

        let mut i = node.first_point();
        let mut j = i + node.num_points - 1;
        while i < j {
            if self.points[i][axis] < split_val {
                i += 1;
            } else {
                self.data.swap(i, j);
                self.points.swap(i, j);
                j -= 1;
            }
        }

        let left_count = i - node.first_point();
        if left_count == 0 || left_count == node.num_points {
            return;
        }

        let node = &mut self.nodes[idx];
        let (old_fst_pt, num_pts) = node.set_left_child(self.nodes_used);
        self.nodes_used += 2;
        let left_child_idx = node.left_child();
        let right_child_idx = node.right_child();

        self.nodes[left_child_idx].set_points(old_fst_pt, left_count);
        self.nodes[right_child_idx].set_points(i, num_pts - left_count);

        self.update_node_bounds(left_child_idx);
        self.update_node_bounds(right_child_idx);

        self.subdivide(left_child_idx, split_kind);
        self.subdivide(right_child_idx, split_kind);
    }
    pub fn nearest(&self, p: &[F; N]) -> (&[F; N], F, ()) {
        const SZ: usize = 32;
        assert_ne!(SZ, 0);
        let mut stack = [&KDNode::EMPTY; SZ];
        let mut stack_ptr = 0;
        macro_rules! push {
            ($n: expr) => {
                assert!(stack_ptr < SZ);
                stack[stack_ptr] = $n;
                stack_ptr += 1;
            };
        }

        macro_rules! pop {
            () => {
                if stack_ptr == 0 {
                    break;
                } else {
                    stack_ptr -= 1;
                    assert_ne!(stack[stack_ptr], &KDNode::EMPTY);
                    stack[stack_ptr]
                }
            };
        }

        let mut curr_best = (0, F::INFINITY);
        let mut node = &self.nodes[self.root_node_idx];
        loop {
            if node.is_leaf() {
                if node.bounds.overlaps(p, curr_best.1).is_none() {
                    continue;
                }
                let fp = node.first_point();
                for i in fp..fp + node.num_points {
                    let pt = self.points[i];
                    let d = dist(&pt, p);
                    if d < curr_best.1 {
                        curr_best = (i, d);
                    }
                }
                node = pop!();
                continue;
            }
            let c1 = &self.nodes[node.left_child()];
            let c2 = &self.nodes[node.right_child()];
            let d1 = c1.bounds.overlaps(p, curr_best.1);
            let d2 = c2.bounds.overlaps(p, curr_best.1);
            node = match (d1, d2) {
                (None, None) => pop!(),
                (None, Some(_)) => c2,
                (Some(_), None) => c1,
                (Some(d1), Some(d2)) => {
                    if d1 < d2 {
                        push!(c2);
                        c1
                    } else {
                        push!(c1);
                        c2
                    }
                }
            };
        }
        let pt = &self.points[curr_best.0];
        let dist = curr_best.1;
        (pt, dist, self.data[curr_best.0])
    }
}

#[test]
fn test_new_kdtree() {
    for n in 1..=100 {
        let pts = (0..n).map(|i| [(i as F).sin(), (i as F).cos()]);
        let kdt = KDTree::<F, (), 2>::new(pts, ());
        let near_to = [0.; 2];
        let found_nearest = kdt.nearest(&near_to).0;

        let pts = (0..n).map(|i| [(i as F).sin(), (i as F).cos()]);
        let naive_nearest = pts
            .map(|p| (p, dist(&near_to, &p)))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .unwrap()
            .0;
        if naive_nearest != *found_nearest {
            let d0 = dist(&near_to, &naive_nearest);
            let d1 = dist(&near_to, &found_nearest);
            assert!((d0 - d1).abs() < 1e-5);
        }
    }
}
