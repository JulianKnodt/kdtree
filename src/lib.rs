#![feature(generic_arg_infer)]
#![feature(test)]
#![feature(extract_if)]
#![feature(let_chains)]
#![feature(trait_alias)]

#[cfg(feature = "usize_idx")]
pub type IdxTy = usize;
#[cfg(not(feature = "usize_idx"))]
pub type IdxTy = u32;

struct Dist<F>(std::marker::PhantomData<F>);

macro_rules! dist_impl {
    ($F: ty) => {
        impl Dist<$F> {
            fn dist_sq<const N: usize>(a: &[$F; N], b: &[$F; N]) -> $F {
                (0..N).map(|i| a[i] - b[i]).map(|v| v * v).sum::<$F>()
            }
            #[inline]
            pub fn dist<const N: usize>(a: &[$F; N], b: &[$F; N]) -> $F {
                Self::dist_sq(a, b).sqrt()
            }
        }
    };
}

dist_impl!(f32);
dist_impl!(f64);

struct SimplexDist<F, const S: usize>(std::marker::PhantomData<F>);

macro_rules! simplex_dist_impl {
    ($F: ty, $math: ident) => {
        impl<const S: usize> SimplexDist<$F, S> {
            fn dist_sq<const N: usize>(a: &[[$F; N]; S], b: &[$F; N]) -> $F {
                match a.as_slice() {
                    [] => 0.,
                    [pt] => Dist::<$F>::dist_sq(&pt, b),
                    // TODO edge?

                    // tri
                    [t0, t1, t2] => {
                        use std::array::from_fn;
                        if N == 3 {
                            $math::tri_sdf_3d(
                                &[t0, t1, t2].map(|v| from_fn(|i| v[i])),
                                from_fn(|i| b[i]),
                            )
                        } else if N == 2 {
                            $math::tri_sdf_2d(
                                &[t0, t1, t2].map(|v| from_fn(|i| v[i])),
                                from_fn(|i| b[i]),
                            )
                        } else {
                            todo!("not sure how to implement this for N != 2 or N != 3");
                        }
                    }
                    _ => todo!(),
                }
            }
            #[inline]
            pub fn dist<const N: usize>(a: &[[$F; N]; S], b: &[$F; N]) -> $F {
                Self::dist_sq(a, b).sqrt()
            }
        }

        mod $math {
            fn sign(x: $F) -> $F {
                if x == 0. {
                    0.
                } else {
                    x.signum()
                }
            }

            fn dot<const N: usize>(a: [$F; N], b: [$F; N]) -> $F {
                (0..N).map(|i| a[i] * b[i]).sum::<$F>()
            }

            fn dot2<const N: usize>(x: [$F; N]) -> $F {
                dot(x, x)
            }

            fn kmul<const N: usize>(k: $F, v: [$F; N]) -> [$F; N] {
                v.map(|v| v * k)
            }

            fn sub<const N: usize>(a: [$F; N], b: [$F; N]) -> [$F; N] {
                std::array::from_fn(|i| a[i] - b[i])
            }

            fn cross([x, y, z]: [$F; 3], [a, b, c]: [$F; 3]) -> [$F; 3] {
                [y * c - z * b, z * a - x * c, x * b - y * a]
            }

            fn cross_2d([x, y]: [$F; 2], [a, b]: [$F; 2]) -> $F {
                x * b - y * a
            }
            pub fn tri_sdf_2d(&[p0, p1, p2]: &[[$F; 2]; 3], p: [$F; 2]) -> $F {
                let e0 = sub(p1, p0);
                let e1 = sub(p2, p1);
                let e2 = sub(p0, p2);
                let v0 = sub(p, p0);
                let v1 = sub(p, p1);
                let v2 = sub(p, p2);
                let [pq0, pq1, pq2] = [[v0, e0], [v1, e1], [v2, e2]]
                    .map(|[v, e]| sub(v, kmul((dot(v, e) / dot2(e)).clamp(0., 1.), e)));
                let s = sign(cross_2d(e0, e2));
                let dx = dot2(pq0).min(dot2(pq1)).min(dot2(pq2));
                let dy = (s * cross_2d(v0, e0))
                    .min(s * cross_2d(v1, e1))
                    .min(s * cross_2d(v2, e2));
                -dx.sqrt() * sign(dy)
            }

            pub fn tri_sdf_3d(&[a, b, c]: &[[$F; 3]; 3], p: [$F; 3]) -> $F {
                let ba = sub(b, a);
                let pa = sub(p, a);
                let cb = sub(c, b);
                let pb = sub(p, b);
                let ac = sub(a, c);
                let pc = sub(p, c);
                let nor = cross(ba, ac);

                let cond = sign(dot(cross(ba, nor), pa))
                    + sign(dot(cross(cb, nor), pb))
                    + sign(dot(cross(ac, nor), pc))
                    < 2.0;
                let v = if cond {
                    let [a, b, c] = [[ba, pa], [cb, pb], [ac, pc]].map(|[e, p]| {
                        let v = kmul((dot(e, p) / dot2(e)).clamp(0., 1.), e);
                        dot2(sub(v, p))
                    });
                    a.min(b).min(c)
                } else {
                    dot(nor, pa) * dot(nor, pa) / dot2(nor)
                };
                v.sqrt()
            }
        }
    };
}

simplex_dist_impl!(f32, math_f32);
simplex_dist_impl!(f64, math_f64);

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
            pub fn sphere(&self) -> Sphere<$F, N> {
                Sphere {
                    center: self.center(),
                    radius: self.extent_length() / 2.,
                }
            }
            pub fn largest_dimension(&self) -> usize {
                (0..N)
                    .map(|i| (i, self.max[i] - self.min[i]))
                    .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
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
            /// If this sphere and another sphere defined by a point and radius overlap, returns
            /// their distance.
            #[inline]
            fn overlaps(&self, pt: &[$F; N], rad: $F) -> Option<$F> {
                let total_d = Dist::<$F>::dist(&self.center, pt);
                debug_assert!(rad >= 0.);
                debug_assert!(self.radius >= 0.);
                let sub_d = total_d - self.radius;
                (sub_d <= rad).then_some(sub_d)
            }

            fn aabb(&self) -> AABB<$F, N> {
                let min = std::array::from_fn(|i| self.center[i] - self.radius);
                let max = std::array::from_fn(|i| self.center[i] + self.radius);
                AABB { min, max }
            }

            fn volume(&self) -> $F {
                const PI: $F = std::f64::consts::PI as $F;
                4. / 3. * PI * self.radius * self.radius * self.radius
            }

            /*
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
            */
        }
    };
}
impl_sphere!(f32);
impl_sphere!(f64);

#[derive(Debug, Copy, Clone, PartialEq)]
struct KDNode<F, const N: usize> {
    bounds: Sphere<F, N>,

    left_child_or_first_point: IdxTy,
    num_points: IdxTy,
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
        self.left_child_or_first_point as usize
    }

    fn is_leaf(&self) -> bool {
        self.num_points > 0
    }

    fn set_left_child(&mut self, left_child: IdxTy) -> (usize, usize) {
        assert!(self.is_leaf());
        let old_first_pt = std::mem::replace(&mut self.left_child_or_first_point, left_child);
        let num_prims = std::mem::take(&mut self.num_points);
        (old_first_pt as usize, num_prims as usize)
    }
    fn first_point(&self) -> usize {
        assert!(self.is_leaf());
        self.left_child_or_first_point as usize
    }
    fn set_points(&mut self, first_pt: IdxTy, num_pts: IdxTy) {
        self.left_child_or_first_point = first_pt;
        self.num_points = num_pts;
    }
}

#[derive(Debug, Clone, Default)]
pub struct KDTree<Q, const N: usize, const S: usize = 1, F = f32> {
    nodes: Vec<KDNode<F, N>>,
    root_node_idx: IdxTy,
    nodes_used: IdxTy,

    elems: Vec<[[F; N]; S]>,
    data: Vec<Q>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SplitKind {
    /// Split along the middle of each sphere
    #[default]
    Midpoint,

    MinMaxVolumeLin(usize),
    // TODO add SAH and other volume heuristic
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum UpdateKind<T> {
    Some(T),
    Delete,
    None,
}

impl<Q, F, const N: usize> KDTree<Q, N, 1, F> {
    pub fn points(&self) -> &[[F; N]] {
        unsafe { std::mem::transmute(self.elems.as_slice()) }
    }
}

macro_rules! impl_kdtree {
    ($F: ty) => {
        impl<const N: usize, const S: usize, T> KDTree<T, N, S, $F> {
            pub fn new(pts: impl Iterator<Item = ([[$F; N]; S], T)>, split: SplitKind) -> Self {
                let (elems, data): (Vec<_>, Vec<_>) = pts.unzip();
                let size = 2 * elems.len() + 1;
                let nodes = vec![KDNode::<$F, N>::EMPTY; size];
                let mut s = Self {
                    nodes,
                    root_node_idx: 0,
                    nodes_used: size.min(1) as IdxTy,
                    elems,
                    data,
                };
                if s.is_empty() {
                    return s;
                }
                s.nodes[0].num_points = s.elems.len() as IdxTy;
                s.update_node_bounds(0);
                s.subdivide(0, split);
                s.nodes.truncate(s.nodes_used as usize);
                s
            }
            fn centroid(&self, i: usize) -> [$F; N] {
                let mut out = [0.; N];
                for pt in self.elems[i] {
                    for i in 0..N {
                        out[i] += pt[i];
                    }
                }
                out.map(|v| v / S as $F)
            }
        }

        impl<const N: usize, const S: usize, T> KDTree<T, N, S, $F> {
            /// Iterate over the data associated with each point
            #[inline]
            pub fn iter_data_mut(&mut self) -> impl Iterator<Item = &mut T> {
                self.data.iter_mut()
            }
            /// Returns the number of points in this KD-tree
            #[inline]
            pub fn len(&self) -> usize {
                self.elems.len()
            }
            /// Returns if there are no points in this KD-tree
            #[inline]
            pub fn is_empty(&self) -> bool {
                self.elems.is_empty()
            }
            fn update_node_bounds(&mut self, idx: usize) {
                let node = &mut self.nodes[idx];
                let mut aabb = AABB::<$F, N>::EMPTY;
                let fp = node.first_point();
                for i in fp..fp + node.num_points as usize {
                    for p in self.elems[i] {
                        aabb.add_point(&p);
                    }
                }
                node.bounds = aabb.sphere();
            }
            fn midpoint_split(&self, node: &KDNode<$F, N>) -> (usize, $F) {
                let mut aabb = AABB::<$F, N>::EMPTY;
                let fp = node.first_point();
                for elem in &self.elems[fp..fp + node.num_points as usize] {
                    for p in elem {
                        aabb.add_point(p);
                    }
                }
                let axis = aabb.largest_dimension();
                (axis, aabb.center()[axis])
            }
            fn min_max_volume_split_lin(&self, node: &KDNode<$F, N>, bins: usize) -> (usize, $F) {
                assert!(bins > 0);
                let mut aabb = AABB::<$F, N>::EMPTY;
                let fp = node.first_point();
                for el in &self.elems[fp..fp + node.num_points as usize] {
                    for p in el {
                        aabb.add_point(&p);
                    }
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
                                for i in fp..fp + node.num_points as usize {
                                    let c = self.centroid(i);
                                    let target_aabb = if c[axis] < split_pt {
                                        &mut left_aabb
                                    } else {
                                        &mut right_aabb
                                    };
                                    for p in self.elems[i] {
                                        target_aabb.add_point(&p);
                                    }
                                }
                                let vol = left_aabb
                                    .sphere()
                                    .volume()
                                    .max(right_aabb.sphere().volume());
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
                        let nl = self.nodes[n.left_child()].bounds.aabb();
                        let nr = self.nodes[n.right_child()].bounds.aabb();
                        self.nodes[i].bounds = nl.add_aabb(&nr).sphere();
                    }
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

                let mut i = node.first_point() as usize;
                let mut j = (i + node.num_points as usize - 1) as usize;
                while i < j {
                    if self.centroid(i)[axis] < split_val {
                        i += 1;
                    } else {
                        self.data.swap(i, j);
                        self.elems.swap(i, j);
                        j -= 1;
                    }
                }

                let left_count = i - node.first_point() as usize;
                if left_count == 0 || left_count == node.num_points as usize {
                    return;
                }

                let node = &mut self.nodes[idx];
                let (old_fst_pt, num_pts) = node.set_left_child(self.nodes_used);
                self.nodes_used += 2;
                let left_child_idx = node.left_child() as usize;
                let right_child_idx = node.right_child() as usize;

                self.nodes[left_child_idx].set_points(old_fst_pt as IdxTy, left_count as IdxTy);
                self.nodes[right_child_idx].set_points(i as IdxTy, (num_pts - left_count) as IdxTy);

                self.update_node_bounds(left_child_idx);
                self.update_node_bounds(right_child_idx);

                self.subdivide(left_child_idx, split_kind);
                self.subdivide(right_child_idx, split_kind);
            }
            /// Returns the nearest point, dist to point, and the associated data with the
            /// point.
            #[inline]
            pub fn nearest(&self, p: &[$F; N]) -> Option<(&[[$F; N]; S], $F, &T)> {
                self.nearest_filter(p, |_| true)
            }
            #[inline]
            pub fn nearest_filter(
                &self,
                p: &[$F; N],
                filter: impl Fn(&T) -> bool,
            ) -> Option<(&[[$F; N]; S], $F, &T)> {
                self.nearest_filter_top_k::<1>(p, <$F>::INFINITY, filter)[0]
            }
            /// Filter allows for skipping elements which return false.
            pub fn nearest_filter_top_k<const K: usize>(
                &self,
                p: &[$F; N],
                ball_radius: $F,
                filter: impl Fn(&T) -> bool,
            ) -> [Option<(&[[$F; N]; S], $F, &T)>; K] {
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
                        unsafe { self.nodes.get_unchecked(n as usize) }
                    }};
                }
                const EMPTY: usize = usize::MAX;

                let mut curr_bests = [(EMPTY, ball_radius); K];
                push!(self.root_node_idx);
                loop {
                    let node = pop!();
                    if node.is_leaf() {
                        let fp = node.first_point();
                        for i in fp..fp + node.num_points as usize {
                            let i = i as usize;
                            if !filter(&self.data[i]) {
                                continue;
                            }
                            let pt = unsafe { self.elems.get_unchecked(i) };
                            let d = SimplexDist::<$F, S>::dist(pt, p);
                            if d < curr_bests[K - 1].1 {
                                curr_bests[K - 1] = (i, d);
                                curr_bests.sort_by(|a, b| a.1.total_cmp(&b.1));
                            }
                            if curr_bests[K - 1].1 == 0. {
                                break;
                            }
                        }
                        continue;
                    }
                    let c1 = unsafe { &self.nodes.get_unchecked(node.left_child() as usize) };
                    let c2 = unsafe { &self.nodes.get_unchecked(node.right_child() as usize) };
                    let d1 = c1.bounds.overlaps(p, curr_bests[K - 1].1);
                    let d2 = c2.bounds.overlaps(p, curr_bests[K - 1].1);

                    match (d1, d2) {
                        (None, None) => {}
                        (None, Some(_)) => push!(node.right_child() as IdxTy),
                        (Some(_), None) => push!(node.left_child() as IdxTy),
                        (Some(d1), Some(d2)) => {
                            if d1 < d2 {
                                push!(node.right_child() as IdxTy);
                                push!(node.left_child() as IdxTy);
                            } else {
                                push!(node.left_child() as IdxTy);
                                push!(node.right_child() as IdxTy);
                            }
                        }
                    };
                }

                curr_bests.map(|(idx, dist)| {
                    (idx != EMPTY).then(|| (&self.elems[idx], dist, &self.data[idx]))
                })
            }
        }
    };
}

impl_kdtree!(f32);
impl_kdtree!(f64);

impl<F: std::fmt::Display + Copy, T, const N: usize> KDTree<T, N, 1, F> {
    pub fn save_as_ply(&self, filename: &str) -> std::io::Result<()> {
        use std::fs::File;
        use std::io::{BufWriter, Write};
        let f = File::create(filename)?;
        let mut f = BufWriter::new(f);
        writeln!(f, "ply")?;
        writeln!(f, "format ascii 1.0")?;
        writeln!(f, "element vertex {}", self.elems.len())?;
        writeln!(f, "property float x")?;
        writeln!(f, "property float y")?;
        writeln!(f, "property float z")?;
        writeln!(f, "end_header")?;
        for &pt in self.points() {
            writeln!(f, "{} {} {}", pt[0], pt[1], pt[2])?;
        }
        Ok(())
    }
}

#[test]
fn test_new_kdtree() {
    let pts = (0..100000)
        .map(|i| [(i as f32).sin(), (i as f32).cos()])
        .map(|p| ([p], ()));
    let kdt = KDTree::<(), 2>::new(pts, Default::default());
    println!("{:?}", kdt.nodes_used);

    let f = kdt.nearest(&[0.1; 2]);
    println!("GOT {f:?}");
}

#[test]
fn test_correct() {
    let probes = [[0.; 2], [0., 1.], [-0.5, -0.5], [1., -0.5]];
    for n in 1..=5000 {
        let pts = (0..n).map(|i| {
            [
                (i as f32 * 131.94 + 4.2451 * n as f32).sin(),
                (i as f32 * 239.73 + 3.18 * n as f32).cos(),
            ]
        });
        let kd_vals = pts.clone().map(|v| ([v], ()));
        let kdt = KDTree::<(), 2>::new(kd_vals, Default::default());
        for probe in probes {
            let (found_nearest, d, _) = kdt.nearest(&probe).unwrap();
            let naive_nearest = pts
                .clone()
                .map(|p| (p, Dist::<f32>::dist(&probe, &p)))
                .min_by(|a, b| a.1.total_cmp(&b.1))
                .unwrap()
                .0;
            assert!((Dist::<f32>::dist(&found_nearest[0], &probe) - d).abs() < 1e-8);
            if naive_nearest != found_nearest[0] {
                let d0 = Dist::<f32>::dist(&naive_nearest, &probe);
                let d1 = Dist::<f32>::dist(&found_nearest[0], &probe);
                assert_eq!(d0, d1);
            }
        }
    }
}

#[test]
fn test_dense_3d() {
    let pts = (0..50000)
        .map(|i| {
            [
                (i as f32 * 131.94 + 4.2451).sin(),
                (i as f32 * 239.73 + 3.18).cos(),
                (i as f32 * 83.38 + 19.32).sin(),
            ]
        })
        .collect::<Vec<_>>();
    let kd_vals = pts.iter().map(|&v| ([v], ()));

    let kdt = KDTree::<(), 3>::new(kd_vals, Default::default());

    let probes = (0..500).map(|i| {
        [
            (i as f32 * 1.1239 + 4.2451).sin(),
            (i as f32 * 2.438 + 3.18).cos(),
            (i as f32 * 8.239 + 19.32).sin(),
        ]
    });
    for p in probes {
        let found_nearest = kdt.nearest(&p).unwrap().0;
        let (naive_nearest, naive_dist) = pts
            .iter()
            .map(|pt| (pt, Dist::<f32>::dist(&p, &pt)))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .unwrap();
        if *naive_nearest != found_nearest[0] {
            assert_eq!(naive_dist, Dist::<f32>::dist(&found_nearest[0], &p));
        }
    }
}

#[cfg(test)]
extern crate test;

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
        ([p], ())
    });

    use core::hint::black_box;
    let kdt = KDTree::<(), _>::new(pts, Default::default());
    let mut i = 0;
    b.iter(|| {
        i += 1;
        let i = i as f32;
        let near_to = [i % 0.3, i % 0.21, i % 0.7, i % 0.893];
        kdt.nearest(&black_box(near_to));
    });
}
