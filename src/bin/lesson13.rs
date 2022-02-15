use ndarray::prelude::*;
use ndarray_npy::write_npy;

fn main() {
    // \pdv[2]{p}{x} + \pdv[2]{p}{y} = b
    let nx = 50;
    let ny = 50;
    let nt = 100;
    let (xmin, xmax) = (0., 2.);
    let (ymin, ymax) = (0., 1.);

    let dx = (xmax - xmin) / (nx as f64 - 1.);
    let dy = (ymax - ymin) / (ny as f64 - 1.);
    let mut p = Array::zeros((ny, nx));
    let mut b = Array::zeros((ny, nx));
    b[[ny / 4, nx / 4]] = 100.;
    b[[3 * ny / 4, 3 * nx / 4]] = -100.;

    for _ in 0..nt {
        let pn = p.clone();
        let b0 = &b.slice(s![1..-1, 1..-1]);
        let p1 = &pn.slice(s![1..-1, 2..]);
        let p2 = &pn.slice(s![1..-1, ..-2]);
        let p3 = &pn.slice(s![2.., 1..-1]);
        let p4 = &pn.slice(s![..-2, 1..-1]);
        let rhs = ((p1 + p2) * dy.powi(2) + (p3 + p4) * dx.powi(2) - b0 * dx.powi(2) * dy.powi(2))
            / (2. * (dx.powi(2) + dy.powi(2)));
        p.slice_mut(s![1..-1, 1..-1]).assign(&rhs);
        p.slice_mut(s![0, ..]).fill(0.);
        p.slice_mut(s![-1, ..]).fill(0.);
        p.slice_mut(s![.., 0]).fill(0.);
        p.slice_mut(s![.., -1]).fill(0.);
    }
    write_npy("array.npy", &p).unwrap();
}
