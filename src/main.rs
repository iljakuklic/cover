
#[cfg(test)]
extern crate quickcheck;
#[cfg(test)]
#[macro_use]
extern crate quickcheck_macros;

mod structs;

/*
use structs::binary_matrix::RowIdx;
use structs::binary_matrix::ColIdx;
use structs::binary_matrix;
*/
//use structs::binary_matrix::{BinaryMatrix, RowIdx, ColIdx};
//use sorted_iter::assume::AssumeSortedByItemExt;
//use std::collections::BTreeSet;

fn sudoku(nb: u32)  {
    let n = nb * nb;
    let item = |x, y, v| {
        format!("R{}C{}V{}", y + 1, x + 1, v + 1)
    };
    for x in 0..n {
        for y in 0..n {
            for v in 0..n {
                print!(" {}", item(x, y, v));
            }
            println!("");
        }
    }
    for x in 0..n {
        for v in 0..n {
            for y in 0..n {
                print!(" {}", item(x, y, v));
            }
            println!("");
        }
    }
    for y in 0..n {
        for v in 0..n {
            for x in 0..n {
                print!(" {}", item(x, y, v));
            }
            println!("");
        }
    }
    for xb in 0..nb {
        for yb in 0..nb {
            for v in 0..n {
                for x in (nb * xb)..(nb * xb + nb) {
                    for y in (nb * yb)..(nb * yb + nb) {
                        print!(" {}", item(x, y, v));
                    }
                }
                println!("");
            }
        }
    }
    //ones.insert(item(col, 1, 3, 2));
    //ones.insert(item(col+1, 2, 2, 2));
    //BinaryMatrix::from_iter(ones.into_iter())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: u32 = args[1].parse().unwrap();
    sudoku(n);
}
