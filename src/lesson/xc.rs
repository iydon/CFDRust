use ndarray::prelude::*;

use crate::macros::{assign, fill, i};

pub fn default() -> (Array2<f64>, Array1<f64>, f64, f64, f64) {
    let nx = 31;
    let ny = 31;

    let nu = 0.05;
    let sigma = 0.25;
    let eps = 1e-4;

    let dx = 2. / (nx as f64 - 1.);
    let dy = 2. / (ny as f64 - 1.);
    let dt = sigma * dx * dy / nu;

    let y = Array::linspace(0., 1., nx);
    let mut p = Array::zeros((ny, nx));
    set_boundary_conditions(&mut p, &y);

    return (p, y, dx, dy, eps);
}

pub fn solve(p: &mut Array2<f64>, y: Array1<f64>, dx: f64, dy: f64, eps: f64) {
    let mut norm = 1.;
    while norm > eps {
        let pn = p.clone();

        assign!(
            p[1..-1, 1..-1] = (dy.powi(2) * (i!(pn[1..-1, 2..]) + i!(pn[1..-1, ..-2]))
                + dx.powi(2) * (i!(pn[2.., 1..-1]) + i!(pn[..-2, 1..-1])))
                / (2. * (dx.powi(2) + dy.powi(2)))
        );

        set_boundary_conditions(p, &y);

        norm = (p.mapv(f64::abs) - pn.mapv(f64::abs)).sum().abs() / pn.mapv(f64::abs).sum();
    }
}

fn set_boundary_conditions(p: &mut Array2<f64>, y: &Array1<f64>) {
    fill!(p[.., 0] = 0.); // p = 0 @ x = 0
    assign!(p[.., -1] = y); // p = y @ x = 2
    let pn = i!(p[1, ..]).to_owned();
    assign!(p[0, ..] = pn); // dp/dy = 0 @ y = 0
    let pn = i!(p[-2, ..]).to_owned();
    assign!(p[-1, ..] = pn); // dp/dy = 0 @ y = 1
}
