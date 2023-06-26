#![feature(generic_arg_infer)]
#![feature(test)]
#![feature(extract_if)]
#![feature(let_chains)]
extern crate test;

use std::collections::{BTreeMap, BTreeSet};

struct Dist<F>(std::marker::PhantomData<F>);

macro_rules! Dist_impl {
    ($F: ty) => {
        impl Dist<$F> {
            fn dist_sq<const N: usize>(a: &[$F; N], b: &[$F; N]) -> $F {
                let mut d = 0.;
                for i in 0..N {
                    let v = a[i] - b[i];
                    d += v * v;
                }
                d
            }
            #[inline]
            pub fn dist<const N: usize>(a: &[$F; N], b: &[$F; N]) -> $F {
                Self::dist_sq(a, b).sqrt()
            }
        }
    };
}

Dist_impl!(f32);
Dist_impl!(f64);

#[derive(Debug, Clone, Copy, PartialEq)]
struct AABB<F, const N: usize> {
    min: [F; N],
    max: [F; N],
}

macro_rules! impl_aabb {
    ($F: ty) => {
        impl<const N: usize> AABB<$F, N> {
            const EMPTY: Self = AABB {
                min: [<$F>::INFINITY; N],
                max: [<$F>::NEG_INFINITY; N],
            };
            pub fn add_point(&mut self, p: &[$F; N]) {
                for i in 0..N {
                    self.min[i] = p[i].min(self.min[i]);
                    self.max[i] = p[i].max(self.max[i]);
                }
            }
            pub fn center(&self) -> [$F; N] {
                std::array::from_fn(|i| (self.min[i] + self.max[i]) / 2.)
            }
            #[inline]
            pub fn extent_length(&self) -> $F {
                Dist::<$F>::dist(&self.max, &self.min)
            }
            #[inline]
            pub fn extent(&self) -> [$F; N] {
                std::array::from_fn(|i| self.max[i] - self.min[i])
            }
            #[inline]
            pub fn to_sphere(&self) -> Sphere<$F, N> {
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
            pub fn add_aabb(&self, o: &Self) -> Self {
                let mut new = *self;
                new.add_point(&o.min);
                new.add_point(&o.max);
                new
            }
        }
    };
}

impl_aabb!(f32);
impl_aabb!(f64);

#[derive(Debug, Clone, Copy, PartialEq)]
struct Sphere<F, const N: usize> {
    center: [F; N],
    radius: F,
}

macro_rules! impl_sphere {
    ($F: ty) => {
        impl<const N: usize> Sphere<$F, N> {
            const EMPTY: Self = Self {
                center: [0.; N],
                radius: <$F>::INFINITY,
            };
            // If it overlaps, returns the distance to the sphere
            #[inline]
            fn overlaps(&self, pt: &[$F; N], rad: $F) -> Option<$F> {
                let d = Dist::<$F>::dist(&self.center, pt);
                debug_assert!(rad >= 0.);
                debug_assert!(self.radius >= 0.);
                let sub_d = d - self.radius;
                (sub_d < rad).then_some(sub_d)
            }

            fn to_aabb(&self) -> AABB<$F, N> {
                let min = std::array::from_fn(|i| self.center[i] - self.radius);
                let max = std::array::from_fn(|i| self.center[i] + self.radius);
                AABB { min, max }
            }

            fn volume(&self) -> $F {
                const PI: $F = std::f64::consts::PI as $F;
                4. / 3. * PI * self.radius * self.radius * self.radius
            }

            #[inline]
            fn contains(&self, p: &[$F; N]) -> bool {
                Dist::<$F>::dist(&self.center, p) < self.radius
            }
            #[inline]
            fn add_point(&mut self, p: &[$F; N]) {
                self.radius = self.radius.max(Dist::<$F>::dist(&self.center, p));
            }
            #[inline]
            fn contains_sphere(&self, s: &Self) -> bool {
                let c_dist = Dist::<$F>::dist(&self.center, &s.center);
                c_dist < (self.radius - s.radius).abs()
            }
            #[inline]
            fn add_sphere(&mut self, s: &Self) {
                let c_dist = Dist::<$F>::dist(&self.center, &s.center);
                self.radius += (if self.radius < s.radius {
                    self.radius - s.radius
                } else {
                    self.radius
                } - c_dist)
                    .max(0.);
            }
        }
    };
}
impl_sphere!(f32);
impl_sphere!(f64);

#[derive(Debug, Copy, Clone, PartialEq)]
struct KDNode<F, const N: usize> {
    bounds: Sphere<F, N>,

    left_child_or_first_point: usize,
    num_points: usize,
}
impl<const N: usize> KDNode<f32, N> {
    const EMPTY: Self = KDNode {
        bounds: Sphere::<f32, N>::EMPTY,
        left_child_or_first_point: 0,
        num_points: 0,
    };
}

impl<const N: usize> KDNode<f64, N> {
    const EMPTY: Self = KDNode {
        bounds: Sphere::<f64, N>::EMPTY,
        left_child_or_first_point: 0,
        num_points: 0,
    };
}

impl<F, const N: usize> KDNode<F, N> {
    #[inline]
    fn right_child(&self) -> usize {
        debug_assert!(!self.is_leaf());
        self.left_child() + 1
    }
    #[inline]
    fn left_child(&self) -> usize {
        debug_assert!(!self.is_leaf());
        self.left_child_or_first_point
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
pub struct KDTree<Q, const N: usize, const ALLOW_UPDATES: bool = false, F = f32> {
    nodes: Vec<KDNode<F, N>>,
    root_node_idx: usize,
    nodes_used: usize,

    points: Vec<[F; N]>,
    data: Vec<Q>,
    // map from data -> (node, index)
    index: BTreeMap<Q, Vec<(usize, usize)>>,

    /// marks which points are no longer valid
    invalids: BTreeSet<usize>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SplitKind {
    /// Split along the middle of each sphere
    Midpoint,

    MinMaxVolumeLin(usize),
    // TODO add SAH and other volume heuristic
}

impl From<()> for SplitKind {
    fn from((): ()) -> SplitKind {
        SplitKind::Midpoint
        //SplitKind::MinMaxVolumeLin(64)
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum UpdateKind<T> {
    Some(T),
    Delete,
    None,
}

macro_rules! impl_kdtree {
    ($F: ty) => {
        impl<const N: usize, T> KDTree<T, N, false, $F> {
            pub fn new(
                pts: impl Iterator<Item = ([$F; N], T)>,
                split: impl Into<SplitKind>,
            ) -> Self {
                let (points, data): (Vec<_>, Vec<_>) = pts.unzip();
                let size = 2 * points.len() + 1;
                let nodes = vec![KDNode::<$F, N>::EMPTY; size];
                let mut s = Self {
                    nodes,
                    root_node_idx: 0,
                    nodes_used: size.min(1),
                    points,
                    data,
                    invalids: BTreeSet::new(),
                    index: BTreeMap::new(),
                };
                if s.is_empty() {
                    return s;
                }
                s.nodes[0].num_points = s.points.len();
                s.update_node_bounds(0);
                s.subdivide(0, split.into());
                s.nodes.truncate(s.nodes_used);
                s
            }
            pub fn updateable(&self) -> KDTree<T, N, true, $F>
            where
                T: Clone,
                $F: Clone,
            {
                KDTree {
                    nodes: self.nodes.clone(),
                    root_node_idx: 0,
                    nodes_used: self.nodes_used,
                    points: self.points.clone(),
                    data: self.data.clone(),
                    invalids: self.invalids.clone(),
                    index: self.index.clone(),
                }
            }
        }
        impl<const N: usize, T, const AU: bool> KDTree<T, N, AU, $F> {
            #[inline]
            pub fn init_index(&mut self)
            where
                T: Ord + Copy,
            {
                for (ni, n) in self.nodes.iter().enumerate() {
                    if n.is_leaf() {
                        let fp = n.first_point();
                        for i in fp..fp + n.num_points {
                            self.index.entry(self.data[i]).or_default().push((ni, i));
                        }
                    }
                }
            }
            /// Rebalances the tree after updates, deleting empty nodes.
            pub fn rebalance(&mut self, new_points: impl Iterator<Item = ([$F; N], T)>)
            where
                T: Ord + Copy,
            {
                if self.is_empty() {
                    return;
                }
                assert!(AU || self.invalids.is_empty());
                while let Some(l) = self.invalids.pop_last() {
                    self.points.swap_remove(l);
                    self.data.swap_remove(l);
                }
                for (p, d) in new_points {
                    self.points.push(p);
                    self.data.push(d);
                }
                assert_eq!(self.points.len(), self.data.len());

                let size = 2 * self.points.len() + 1;
                for n in &mut self.nodes {
                    *n = KDNode::<$F, N>::EMPTY;
                }

                self.nodes.resize(size, KDNode::<$F, N>::EMPTY);

                self.nodes_used = 1;
                self.nodes[0].num_points = self.points.len();
                self.update_node_bounds(0);
                self.subdivide(0, ().into());
                self.nodes.truncate(self.nodes_used);
                if !self.index.is_empty() {
                    self.index = BTreeMap::new();
                    self.init_index();
                }
            }
            /// Iterate over the data associated with each point
            #[inline]
            pub fn iter_data_mut(&mut self) -> impl Iterator<Item = &mut T> {
                self.data.iter_mut()
            }
            /// Returns the number of points in this KD-tree
            #[inline]
            pub fn len(&self) -> usize {
                self.points.len()
            }
            /// Returns if there are no points in this KD-tree
            #[inline]
            pub fn is_empty(&self) -> bool {
                self.points.is_empty()
            }
            fn update_node_bounds(&mut self, idx: usize) {
                let node = &mut self.nodes[idx];
                let mut aabb = AABB::<$F, N>::EMPTY;
                let fp = node.first_point();
                for i in fp..fp + node.num_points {
                    let p = &self.points[i];
                    if AU && self.invalids.contains(&i) {
                        continue;
                    }
                    aabb.add_point(p);
                }
                node.bounds = aabb.to_sphere();
            }
            fn midpoint_split(&self, node: &KDNode<$F, N>) -> (usize, $F) {
                let mut aabb = AABB::<$F, N>::EMPTY;
                let fp = node.first_point();
                for p in &self.points[fp..fp + node.num_points] {
                    aabb.add_point(p);
                }
                let axis = aabb.largest_dimension();
                (axis, aabb.center()[axis])
            }
            fn min_max_volume_split_lin(&self, node: &KDNode<$F, N>, bins: usize) -> (usize, $F) {
                assert!(bins > 0);
                let mut aabb = AABB::<$F, N>::EMPTY;
                let fp = node.first_point();
                for p in &self.points[fp..fp + node.num_points] {
                    aabb.add_point(p);
                }
                let (axis, best_pos, _) = (0..N)
                    .map(|axis| {
                        let (best_pos, min_max_vol) = (0..bins)
                            .map(|i| {
                                let frac = (i as $F) / (bins as $F);
                                let split_pt = aabb.min[axis] + frac * aabb.extent()[axis];
                                let mut left_aabb = AABB::<$F, N>::EMPTY;
                                let mut right_aabb = AABB::<$F, N>::EMPTY;
                                let fp = node.first_point();
                                for p in &self.points[fp..fp + node.num_points] {
                                    if p[axis] < split_pt {
                                        left_aabb.add_point(p);
                                    } else {
                                        right_aabb.add_point(p);
                                    }
                                }
                                let vol = left_aabb
                                    .to_sphere()
                                    .volume()
                                    .max(right_aabb.to_sphere().volume());
                                (split_pt, vol)
                            })
                            .min_by(|a, b| a.1.total_cmp(&b.1))
                            .unwrap();
                        (axis, best_pos, min_max_vol)
                    })
                    .min_by(|a, b| a.2.total_cmp(&b.2))
                    .unwrap();

                (axis, best_pos)
            }
            pub fn refit(&mut self) {
                // note do not need to do multiple iterations
                // because the children are always greater than the parents in index
                for i in (0..self.nodes.len()).rev() {
                    let n = &self.nodes[i];
                    if n.is_leaf() {
                        self.update_node_bounds(i)
                    } else {
                        let nl = self.nodes[n.left_child()].bounds.to_aabb();
                        let nr = self.nodes[n.right_child()].bounds.to_aabb();
                        self.nodes[i].bounds = nl.add_aabb(&nr).to_sphere();
                    }
                }
            }
            pub fn adjust_points_matching(
                &mut self,
                matching: &[T],
                mut adj: impl FnMut([$F; N], T) -> UpdateKind<([$F; N], T)>,
            ) where
                T: Copy + Ord,
            {
                assert!(!self.index.is_empty());

                for &m in matching {
                    let Some(mut idxs) = self.index.remove(&m) else {
                        continue;
                    };
                    let iter =
                        idxs.extract_if(|&mut (ni, i)| match adj(self.points[i], self.data[i]) {
                            UpdateKind::Some((p, d)) => {
                                // TODO maybe eagerly update parents here?
                                if !self.nodes[ni].bounds.contains(&p) {
                                    self.nodes[ni].bounds.add_point(&p);
                                    let mut parent = ni.checked_sub(2);
                                    // TODO
                                    while let Some(p) = parent &&
                                !self.nodes[p].bounds.contains_sphere(&self.nodes[p+2].bounds) {
                                    let s = self.nodes[p+2].bounds;
                                    self.nodes[p].bounds.add_sphere(&s);
                                    parent = p.checked_sub(2);
                                }
                                }
                                assert!(!self.invalids.contains(&i));
                                self.points[i] = p;
                                self.data[i] = d;
                                if d != m {
                                    self.index.entry(d).or_default().push((ni, i));
                                    true
                                } else {
                                    false
                                }
                            }
                            UpdateKind::Delete => {
                                self.invalids.insert(i);
                                true
                            }
                            UpdateKind::None => false,
                        });

                    // drain iter to ensure it's actually filtered.
                    for _ in iter {}

                    assert_eq!(self.index.insert(m, idxs), None);
                }
            }
            fn subdivide(&mut self, idx: usize, split_kind: SplitKind) {
                let node = &self.nodes[idx];
                // TODO here can use a different amount of points so it will be faster?
                if node.num_points <= 8 {
                    return;
                }
                let (axis, split_val) = match split_kind {
                    SplitKind::Midpoint => self.midpoint_split(node),
                    SplitKind::MinMaxVolumeLin(bins) => self.min_max_volume_split_lin(node, bins),
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
            #[inline]
            pub fn nearest(&self, p: &[$F; N]) -> Option<(&[$F; N], $F, &T)> {
                self.nearest_filter(p, |_| true)
            }
            #[inline]
            pub fn nearest_filter(
                &self,
                p: &[$F; N],
                filter: impl Fn(&T) -> bool,
            ) -> Option<(&[$F; N], $F, &T)> {
                self.nearest_filter_top_k::<1>(p, <$F>::INFINITY, filter)[0]
            }
            /// Filter allows for skipping elements which return false.
            pub fn nearest_filter_top_k<const K: usize>(
                &self,
                p: &[$F; N],
                ball_radius: $F,
                filter: impl Fn(&T) -> bool,
            ) -> [Option<(&[$F; N], $F, &T)>; K] {
                if K == 0 {
                    // Need K for type checking
                    return [None; K];
                }
                if self.is_empty() {
                    return [None; K];
                }
                let mut heap = vec![];
                const SZ: usize = 32;
                let mut stack = [0; SZ];
                let mut stack_ptr = 0;
                macro_rules! push {
                    ($n: expr) => {{
                        if stack_ptr == SZ {
                            heap.push($n);
                        } else {
                            unsafe {
                                *stack.get_unchecked_mut(stack_ptr) = $n;
                            }
                            stack_ptr += 1;
                        }
                    }};
                }

                macro_rules! pop {
                    () => {{
                        let n = if let Some(n) = heap.pop() {
                            n
                        } else if stack_ptr == 0 {
                            break;
                        } else {
                            stack_ptr -= 1;
                            unsafe { *stack.get_unchecked(stack_ptr) }
                        };
                        unsafe { self.nodes.get_unchecked(n) }
                    }};
                }
                const EMPTY: usize = usize::MAX;

                let mut curr_bests = [(EMPTY, ball_radius); K];
                push!(self.root_node_idx);
                loop {
                    let node = pop!();
                    if node.is_leaf() {
                        let fp = node.first_point();
                        for i in fp..fp + node.num_points {
                            if AU && self.invalids.contains(&i) {
                                continue;
                            }
                            if !filter(&self.data[i]) {
                                continue;
                            }
                            let pt = unsafe { self.points.get_unchecked(i) };
                            let d = Dist::<$F>::dist(pt, p);
                            if d < curr_bests[K - 1].1 {
                                curr_bests[K - 1] = (i, (d - 1e-9).max(0.));
                                curr_bests.sort_by(|a, b| a.1.total_cmp(&b.1));
                            }
                            if curr_bests[K - 1].1 == 0. {
                                break;
                            }
                        }
                        continue;
                    }
                    let c1 = unsafe { &self.nodes.get_unchecked(node.left_child()) };
                    let c2 = unsafe { &self.nodes.get_unchecked(node.right_child()) };
                    let d1 = c1.bounds.overlaps(p, curr_bests[K - 1].1);
                    let d2 = c2.bounds.overlaps(p, curr_bests[K - 1].1);

                    match (d1, d2) {
                        (None, None) => {}
                        (None, Some(_)) => push!(node.right_child()),
                        (Some(_), None) => push!(node.left_child()),
                        (Some(d1), Some(d2)) => {
                            if d1 < d2 {
                                push!(node.right_child());
                                push!(node.left_child());
                            } else {
                                push!(node.left_child());
                                push!(node.right_child());
                            }
                        }
                    };
                }

                curr_bests.map(|(idx, dist)| {
                    (idx != EMPTY).then(|| (&self.points[idx], dist, &self.data[idx]))
                })
            }
        }
    };
}

impl_kdtree!(f32);
impl_kdtree!(f64);

impl<F: std::fmt::Display + Copy, T, const N: usize, const AU: bool> KDTree<T, N, AU, F> {
    pub fn save_as_ply(&self, filename: &str) -> std::io::Result<()> {
        use std::fs::File;
        use std::io::{BufWriter, Write};
        let f = File::create(filename)?;
        let mut f = BufWriter::new(f);
        writeln!(f, "ply")?;
        writeln!(f, "format ascii 1.0")?;
        writeln!(f, "element vertex {}", self.points.len())?;
        writeln!(f, "property float x")?;
        writeln!(f, "property float y")?;
        writeln!(f, "property float z")?;
        writeln!(f, "end_header")?;
        for &pt in &self.points {
            writeln!(f, "{} {} {}", pt[0], pt[1], pt[2])?;
        }
        Ok(())
    }
    #[inline]
    pub fn points_iter(&self) -> impl Iterator<Item = &[F; N]> {
        self.points.iter()
    }
}

#[test]
fn test_new_kdtree() {
    let pts = (0..100000)
        .map(|i| [(i as f32).sin(), (i as f32).cos()])
        .map(|p| (p, ()));
    let kdt = KDTree::<(), 2>::new(pts, ());
    println!("{:?}", kdt.nodes_used);

    let f = kdt.nearest(&[0.1; 2]);
    println!("GOT {f:?}");
}

#[test]
fn test_correct() {
    for n in 1..=100 {
        let pts = (0..n).map(|i| ([(i as f32).sin(), (i as f32).cos()], ()));
        let kdt = KDTree::<(), 2>::new(pts, ());
        let near_to = [0.; 2];
        let found_nearest = kdt.nearest(&near_to).unwrap().0;

        let pts = (0..n).map(|i| [(i as f32).sin(), (i as f32).cos()]);
        let naive_nearest = pts
            .map(|p| (p, Dist::<f32>::dist(&near_to, &p)))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .unwrap()
            .0;
        if naive_nearest != *found_nearest {
            let d0 = Dist::<f32>::dist(&near_to, &naive_nearest);
            let d1 = Dist::<f32>::dist(&near_to, &found_nearest);
            assert!((d0 - d1).abs() < 1e-5);
        }
    }
}

#[bench]
fn bench_kdtree(b: &mut test::Bencher) {
    let n = 10000000;
    let pts = (0..n).map(|i| i as f32).map(|i| {
        let p = [
            i.sin(),
            (i * 0.07).cos(),
            (i * 0.3).sin(),
            (i * 0.12).cos(), /**/
        ];
        (p, ())
    });

    use core::hint::black_box;
    let kdt = KDTree::<(), _>::new(pts, ());
    let mut i = 0;
    b.iter(|| {
        i += 1;
        let i = i as f32;
        let near_to = [i % 0.3, i % 0.21, i % 0.7, i % 0.893];
        kdt.nearest(&black_box(near_to));
    });
}
