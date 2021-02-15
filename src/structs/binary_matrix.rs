//! A space-efficient implementation of a sparse binary matrix.

use phantom_newtype::Id;
use sorted_iter::{SortedIterator, sorted_iterator::SortedByItem};
use sorted_iter::assume::*;

/// Type for matrix row/column indices.
type Idx<I> = Id<I, u32>;

/// Row index tag
#[derive(Debug, PartialEq, Eq)]
pub enum RowTag {}

/// Column index tag
#[derive(Debug, PartialEq, Eq)]
pub enum ColTag {}

/// Row index, differentiated at type level.
pub type RowIdx = Idx<RowTag>;
/// Column index, differentiated at type level.
pub type ColIdx = Idx<ColTag>;

/// General SCS/CCS sparse binary matric data storage.
///
/// Parametrized by major/minor dimension index tag pantom types
/// for greater type safety.
#[derive(Debug, PartialEq, Eq)]
struct Store<J, I> {
    /// Where the major dimension breaks occur.
    breaks: Vec<usize>,
    /// Minor dimension indices with ones in them.
    ones: Vec<Idx<I>>,
    /// Something to use up the major index marker tag (phantom type).
    major_tag: Id<J, ()>,
}

impl<J, I> Store<J, I> {

    /// New empty matrix data store.
    fn new_empty() -> Store<J, I> {
        Self::from_iter(std::iter::empty::<(Idx<J>, Idx<I>)>())
    }

    /// Create a new matrix store using a sorted iterator over one-cell
    /// coordinates.
    fn from_iter<Iter>(iter: Iter) -> Store<J, I>
        where Iter: SortedIterator<Item = (Idx<J>, Idx<I>)> {
        let mut breaks = Vec::new();
        let mut ones = Vec::new();

        for (b, i) in iter {
            while breaks.len() <= *b.get() as usize {
                breaks.push(ones.len());
            }
            ones.push(i);
        }
        breaks.push(ones.len());
        Self { breaks, ones, major_tag: Id::new(()), }
    }

    /// Size of the store in major dimension.
    fn major_size(&self) -> u32 {
        (self.breaks.len() - 1) as u32
    }

    /// Size of the store in minor dimension. Slow.
    fn minor_size(&self) -> u32 {
        self.span_iter().map(|s| s.last().map_or(0, |&x| x.get() + 1))
            .max().unwrap_or(0)
    }

    /// Iterator over indices in major dimension,
    fn major_idxs(&self) -> impl Iterator<Item = Idx<J>> {
        (0..(self.breaks.len() - 1) as u32).map(Id::new)
    }

    /// Get vector at specified major dimension.
    fn get_major<'a>(&'a self, idx: Idx<J>) -> &'a[Idx<I>] {
        let j = *idx.get() as usize;
        &self.ones[self.breaks[j] as usize .. self.breaks[j+1] as usize]
    }

    /// Iterate over vectors (represented as slices) that constitute the matrix.
    fn span_iter<'a>(&'a self) -> impl Iterator<Item = &'a[Idx<I>]> + 'a {
        self.major_idxs().map(move |j| self.get_major(j))
    }

    /// Iterator over one-cells.
    fn ones_iter<'a>(&'a self)
            -> impl Iterator<Item = (Idx<J>, Idx<I>)> + SortedByItem + 'a {
        self.major_idxs().flat_map(
            move |j| self.get_major(j).iter().map(move |i| (j, *i))
        ).assume_sorted_by_item()
    }

    /// Iterate over one-cells in minor order, effectively transposing the matrix.
    fn transposed_cells_iter<'a>(&'a self)
            -> impl SortedIterator<Item = (Idx<I>, Idx<J>)> + 'a {
        sorted_iter::multiway_union(self.major_idxs().map(move |j| {
            self.get_major(j).iter().map(move |i| (*i, j)).assume_sorted_by_item()
        }))
    }

    /// Create a transposed matrix.
    fn transposed(&self) -> Store<I, J> {
        Store::from_iter(self.transposed_cells_iter())
    }

    const BRAILLE_START: char = '\u{2800}';

    /// Use braille unicode characters to dump the matrix into a String.
    fn dump(&self) -> String {
        let chcols = (self.minor_size() + 1) / 2;
        let empty: &[Idx<I>] = &[];
        // rows with true/false iterators
        let mut rows = self.span_iter().chain(std::iter::repeat(empty)).map(|s| {
            let mut ones = s.iter().copied().peekable();
            let mut iter = (0u32..).map(move |i| {
                let hit = ones.peek() == Some(&Idx::new(i));
                if hit { ones.next(); true } else { false }
            });
            // Group by two bools.
            (0..chcols).map(move |_| (iter.next().unwrap(), iter.next().unwrap()))
        });
        // groups of 4 rows of true/false iterators
        let rows4 = (0..((self.major_size() + 3) / 4)).map(move |_| {
            let r0 = rows.next().unwrap();
            let r1 = rows.next().unwrap();
            let r2 = rows.next().unwrap();
            let r3 = rows.next().unwrap();
            r0.zip(r1).zip(r2).zip(r3).map(|(((a, b), c), d)| (a, b, c, d))
        });
        // Convert the 4x2 entries into chars
        rows4.flat_map(|r4| r4.map(|((x0, x3), (x1, x4), (x2, x5), (x6, x7))| {
            let bits = [x7, x6, x5, x4, x3, x2, x1, x0];
            let offset = bits.iter().fold(0, |x, &d| 2 * x + d as u32);
            std::char::from_u32(Self::BRAILLE_START as u32 + offset).unwrap()
        }).chain(std::iter::once('\n'))).collect()
    }

    /// Drop certain one-cells according to given predicate.
    fn filtered<Pred: Fn(&(Idx<J>, Idx<I>)) -> bool>(&self, p: Pred) -> Self {
        Self::from_iter(self.ones_iter().filter(p))
    }

    /// Check invariants.
    #[cfg(test)]
    fn valid(&self) -> Result<(), &'static str> {
        fn pairwise<It, F>(i: It, f: F) -> bool
        where It: Iterator + Clone, It::Item: Ord,
              F: Fn(It::Item, It::Item) -> bool {
            i.clone().zip(i.skip(1)).all(|(a, b)| f(a, b))
        }
        fn check(b: bool, err: &'static str) -> Result<(), &str> {
            if b { Ok(()) } else { Err(err) }
        }
        check(!self.breaks.is_empty(), "break list empty")?;
        check(self.breaks[0] == 0, "does not start at 0")?;
        check(*self.breaks.last().unwrap() as usize == self.ones.len(),
              "the last index does not match")?;
        check(self.span_iter().map(|s| s.len()).sum::<usize>() == self.ones.len(),
              "span lengths do not match with total size")?;
        check(pairwise(self.breaks.iter(), |a, b| a <= b), "breaks not sorted")?;
        check(self.span_iter().all(|s| pairwise(s.iter(), |a, b| a < b)),
              "spans not sorted")?;
        Ok(())
    }
}

impl<J, I> Clone for Store<J, I> {
    fn clone(&self) -> Self {
        Store {
            breaks: self.breaks.clone(),
            ones: self.ones.clone(),
            major_tag: Id::new(()),
        }
    }
}

/// Sparse binary matrix with efficient iteration over one-cells in both
/// rows and columns.
///
/// The efficient iteration is facilitated by storing matrix both by rows
/// and by colums, i.e. its direct and transposed variants are both stored.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct BinaryMatrix {
    by_row: Store<RowTag, ColTag>,
    by_col: Store<ColTag, RowTag>,
}

impl BinaryMatrix {
    /// Create a binary matrix from a sorted iterator over one-cell coordinates.
    pub fn from_iter<Iter>(it: Iter) -> Self 
    where Iter: SortedIterator<Item = (RowIdx, ColIdx)> {
        Self::from_rows(Store::from_iter(it))
    }

    /// Create a binary matrix from rows store.
    fn from_rows(by_row: Store<RowTag, ColTag>) -> Self {
        let by_col = by_row.transposed();
        Self { by_row, by_col }
    }

    /*
    pub fn cols(self: &Self) -> impl Iterator<Item = &[RowIdx]> {
        self.by_col.span_iter()
    }
    */

}

/// Display BinaryMatrix using Braille patterns.
impl std::fmt::Display for BinaryMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.by_row.dump())
    }
}

#[cfg(test)]
mod test {

    use super::{Store, BinaryMatrix, RowTag, ColTag, RowIdx, ColIdx, Idx};
    use quickcheck::{Arbitrary, Gen};
    use sorted_iter::assume::AssumeSortedByItemExt;
    use rand::Rng;
    use std::boxed::Box;
    use std::hash::{Hash, Hasher};
    use std::collections::hash_map::DefaultHasher;
    use std::collections::btree_set::BTreeSet;

    // Arbitrary instances

    impl<J: 'static, I: 'static> Arbitrary for Store<J, I>
    where J: std::marker::Send, I: std::marker::Send {
        fn arbitrary<G: Gen>(g: &mut G) -> Self {
            let nmaj = g.gen_range(0u32, g.size() as u32);
            let nmin = g.gen_range(0u32, g.size() as u32);
            let density = g.gen::<f32>().powi(4);
            let entries = (0..nmaj).flat_map(|j| (0..nmin).map(
                    move |i| (Idx::new(j), Idx::new(i))))
                .filter(|_| g.gen::<f32>() < density);
            Store::from_iter(entries.assume_sorted_by_item())
        }

        fn shrink(&self) -> Box<dyn Iterator<Item = Self>> {
            let orig = self.clone();
            let shrink_major = self.major_size().shrink()
                .map(move |n| orig.clone().filtered(|&(j, _)| *j.get() < n));
            let orig = self.clone();
            let shrink_minor = self.minor_size().shrink()
                .map(move |n| orig.filtered(|&(_, i)| *i.get() < n));
            let cells: BTreeSet<_> = self.ones_iter()
                .map(|(j, i)| (*j.get(), *i.get())).collect();
            let shrink_cells = cells.shrink().map(|cs| Self::from_iter(
                    cs.iter().map(|(j, i)| (Idx::new(*j), Idx::new(*i)))
                        .assume_sorted_by_item()));
            if self.major_size() >= self.minor_size() {
                Box::new(shrink_major.chain(shrink_minor).chain(shrink_cells))
            } else {
                Box::new(shrink_minor.chain(shrink_major).chain(shrink_cells))
            }
        }
    }

    impl Arbitrary for BinaryMatrix {
        fn arbitrary<G: Gen>(g: &mut G) -> Self {
            Self::from_rows(Store::arbitrary(g))
        }

        fn shrink(&self) -> Box<dyn Iterator<Item = Self>> {
            Box::new(self.by_row.shrink().map(Self::from_rows))
        }
    }

    // Helper functions

    fn hash_predicate<T: Hash>(x: &T) -> bool {
        let mut h = DefaultHasher::new();
        x.hash(&mut h);
        h.finish() & 4 == 1
    }

    fn make_store_unchecked<I>(i: I) -> Store<RowTag, ColTag>
    where I: Iterator<Item = (u32, u32)> {
        Store::from_iter(i.map(|(r, c)| (RowIdx::new(r), ColIdx::new(c)))
                          .assume_sorted_by_item())
    }

    // Unit tests for corner cases

    #[test]
    fn unit_empty_valid() -> Result<(), &'static str> {
        Store::<RowIdx, ColIdx>::new_empty().valid()
    }

    #[test]
    fn unit_width0() {
        assert_eq!(Store::<RowIdx, ColIdx>::new_empty().minor_size(), 0);
    }

    #[test]
    fn unit_width2() {
        let s1 = [(0, 0), (3, 1)].iter().copied();
        assert_eq!(make_store_unchecked(s1).minor_size(), 2);
    }

    #[test]
    fn unit_dump() {
        let entries = [
            (vec![], ""),
            (vec![(0, 0)], "⠁\n"),
            (vec![(0, 1)], "⠈\n"),
            (vec![(1, 0)], "⠂\n"),
            (vec![(1, 1)], "⠐\n"),
            (vec![(2, 0)], "⠄\n"),
            (vec![(2, 1)], "⠠\n"),
            (vec![(3, 0)], "⡀\n"),
            (vec![(3, 1)], "⢀\n"),
            (vec![(1, 3)], "⠀⠐\n"),
            (vec![(0, 0), (0, 1)], "⠉\n"),
            (vec![(0, 0), (1, 0)], "⠃\n"),
            (vec![(0, 0), (1, 1)], "⠑\n"),
            (vec![(0, 0), (3, 1)], "⢁\n"),
            (vec![(0, 7), (1, 3)], "⠀⠐⠀⠈\n"),
            (vec![(3, 1), (7, 0)], "⢀\n⡀\n"),
        ];
        for (s, d) in entries.iter() {
            let s = make_store_unchecked(s.iter().copied());
            assert_eq!(s.dump(), *d);
        }
    }

    // Test that generator and shrinker produce valid matrices

    #[quickcheck]
    fn qc_arbitrary_valid(s: Store<RowIdx, ColIdx>) -> Result<(), &'static str> {
        s.valid()
    }

    #[quickcheck]
    fn qc_shrinked_valid(s: Store<RowIdx, ColIdx>) -> Result<(), &'static str> {
        for ss in s.shrink().take(30) { ss.valid()?; }
        Ok(())
    }

    // Check operations preserve invariants

    #[quickcheck]
    fn qc_filtered_valid(s: Store<RowIdx, ColIdx>) -> Result<(), &'static str> {
        s.filtered(hash_predicate).valid()
    }

    #[quickcheck]
    fn qc_transposed_valid(s: Store<RowIdx, ColIdx>) -> Result<(), &'static str> {
        s.transposed().valid()
    }

    // Check some basic properties

    #[quickcheck]
    fn qc_filtered_trivial_true(s: Store<RowIdx, ColIdx>) -> bool {
        s.filtered(|_| true) == s
    }

    #[quickcheck]
    fn qc_filtered_trivial_false(s: Store<RowIdx, ColIdx>) -> bool {
        s.filtered(|_| false).major_size() == 0
    }

    #[quickcheck]
    fn qc_filtered_model(s: Store<RowIdx, ColIdx>) -> bool {
        let model: BTreeSet<_> = s.ones_iter().filter(hash_predicate).collect();
        let result: BTreeSet<_> = s.filtered(hash_predicate).ones_iter().collect();
        result == model
    }

    #[quickcheck]
    fn qc_transposed_twice(s: Store<RowIdx, ColIdx>) -> bool {
        s.transposed().transposed() == s
    }

    #[quickcheck]
    fn qc_transposed_model(s: Store<RowIdx, ColIdx>) -> bool {
        let model: BTreeSet<_> = s.ones_iter().map(|(a,b)| (b,a)).collect();
        let result: BTreeSet<_> = s.transposed().ones_iter().collect();
        result == model
    }

    #[quickcheck]
    fn qc_dump_popcount(s: Store<RowIdx, ColIdx>) -> bool {
        let bs = Store::<RowIdx, ColIdx>::BRAILLE_START as u32;
        s.dump().chars().filter_map(|c| {
            let c = c as u32;
            if c & !0xFFu32 == bs {
                Some((c - bs).count_ones())
            } else {
                None
            }
        }).sum::<u32>() == s.ones.len() as u32
    }

    // Tests for BianryMatrix

    #[quickcheck]
    fn qc_bycol_transposed(bm: BinaryMatrix) -> bool {
        bm.by_col.transposed() == bm.by_row
    }

}
