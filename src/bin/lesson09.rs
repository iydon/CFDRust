use ndarray::prelude::*;
use ndarray_npy::write_npy;

fn main() {
    // \pdv{u}{t} = ν\pdv[2]{u}{x} + v\pdv[2]{u}{y}
    let nx = 31;
    let ny = 31;
    let nt = 17;
    let nu = 0.05;
    let sigma = 0.25;

    let dx = 2. / (nx as f64 - 1.);
    let dy = 2. / (ny as f64 - 1.);
    let dt = sigma * dx * dy / nu;
    let mut u = Array::ones((ny, nx));
    let begin = (0.5 / dx) as usize;
    let end = (1. / dx + 1.) as usize;
    u.slice_mut(s![begin..end, begin..end]).fill(2.);  // nx == ny

    for _ in 0..nt + 1 {
        let un = u.clone();
        let x0 = &un.slice(s![1..-1, 1..-1]);
        let x1 = &un.slice(s![1..-1, 2..]);
        let x2 = &un.slice(s![1..-1, ..-2]);
        let x3 = &un.slice(s![2.., 1..-1]);
        let x4 = &un.slice(s![..-2, 1..-1]);
        let rhs = x0
            + nu * dt / dx.powi(2) * (x1 - 2. * x0 + x2)
            + nu * dt / dy.powi(2) * (x3 - 2. * x0 + x4);
        u.slice_mut(s![1..-1, 1..-1]).assign(&rhs);
        u.slice_mut(s![0, ..]).fill(1.);
        u.slice_mut(s![-1, ..]).fill(1.);
        u.slice_mut(s![.., 0]).fill(1.);
        u.slice_mut(s![.., -1]).fill(1.);
    }
    write_npy("array.npy", &u).unwrap();
}
