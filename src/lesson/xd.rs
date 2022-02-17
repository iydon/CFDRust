use ndarray::prelude::*;

use crate::macros::{assign, fill, i};

pub fn default() -> (Array2<f64>, Array2<f64>, usize, f64, f64) {
    let nx = 50;
    let ny = 50;
    let nt = 100;

    let dx = 2. / (nx as f64 - 1.);
    let dy = 1. / (ny as f64 - 1.);

    let mut p = Array::zeros((ny, nx));
    let mut b = Array::zeros((ny, nx));
    b[[ny / 4, nx / 4]] = 100.;
    b[[3 * ny / 4, 3 * nx / 4]] = -100.;

    return (p, b, nt, dx, dy);
}

pub fn solve(p: &mut Array2<f64>, b: &Array2<f64>, nt: usize, dx: f64, dy: f64) {
    for _ in 0..nt + 1 {
        let pn = p.clone();

        assign!(
            p[1..-1, 1..-1] = ((i!(pn[1..-1, 2..]) + i!(pn[1..-1, ..-2])) * dy.powi(2)
                + (i!(pn[2.., 1..-1]) + i!(pn[..-2, 1..-1])) * dx.powi(2)
                - i!(b[1..-1, 1..-1]) * dx.powi(2) * dy.powi(2))
                / (2. * (dx.powi(2) + dy.powi(2)))
        );

        fill!(p[0, ..] = 0.);
        fill!(p[-1, ..] = 0.);
        fill!(p[.., 0] = 0.);
        fill!(p[.., -1] = 0.);
    }
}
