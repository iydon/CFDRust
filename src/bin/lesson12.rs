use ndarray::prelude::*;
use ndarray_npy::write_npy;

fn set_boundary_conditions(p: &mut Array2<f64>, y: &Array1<f64>) {
    p.slice_mut(s![.., 0]).fill(0.); // p = 0 @ x = 0
    p.slice_mut(s![.., -1]).assign(&y); // p = y @ x = 2
    let p1 = p.slice(s![1, ..]).to_owned();
    p.slice_mut(s![0, ..]).assign(&p1); // dp/dy = 0 @ y = 0
    let p2 = p.slice(s![-2, ..]).to_owned();
    p.slice_mut(s![-1, ..]).assign(&p2); // dp/dy = 0 @ y = 1
}

fn laplace2d(p: &mut Array2<f64>, y: &Array1<f64>, dx: f64, dy: f64, eps: f64) {
    let mut l1norm = 1.;
    while l1norm > eps {
        let pn = p.clone();
        let p1 = &pn.slice(s![1..-1, 2..]);
        let p2 = &pn.slice(s![1..-1, ..-2]);
        let p3 = &pn.slice(s![2.., 1..-1]);
        let p4 = &pn.slice(s![..-2, 1..-1]);
        let rhs =
            (dy.powi(2) * (p1 + p2) + dx.powi(2) * (p3 + p4)) / (2. * (dx.powi(2) + dy.powi(2)));
        p.slice_mut(s![1..-1, 1..-1]).assign(&rhs);

        set_boundary_conditions(p, &y);

        l1norm = (p.mapv(f64::abs) - pn.mapv(f64::abs)).sum().abs() / pn.mapv(f64::abs).sum();
    }
}

fn main() {
    // \pdv[2]{p}{x} + \pdv[2]{p}{y} = 0
    let nx = 31;
    let ny = 31;

    let dx = 2. / (nx as f64 - 1.);
    let dy = 1. / (ny as f64 - 1.);
    let y = Array::linspace(0., 1., nx);
    let mut p = Array::zeros((ny, nx));
    set_boundary_conditions(&mut p, &y);
    laplace2d(&mut p, &y, dx, dy, 1e-4);
    write_npy("array.npy", &p).unwrap();
}
