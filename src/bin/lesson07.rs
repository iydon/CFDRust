use ndarray::prelude::*;
use ndarray_npy::write_npy;

fn main() {
    // \pdv{u}{t} + c\pdv{u}{x} + c\pdv{u}{y} = 0
    let nx = 81;
    let ny = 81;
    let nt = 100;
    let c = 1.;
    let sigma = 0.2;

    let dx = 2. / (nx as f64 - 1.);
    let dy = 2. / (ny as f64 - 1.);
    let dt = sigma * dx;
    let mut u = Array::ones((nx, ny));
    let begin = (0.5 / dx) as usize;
    let end = (1. / dx + 1.) as usize;
    u.slice_mut(s![begin..end, begin..end]).fill(2.);

    for _ in 0..nt + 1 {
        let un = u.clone();
        let x = &un.slice(s![1.., 1..]);
        let y = &un.slice(s![1.., ..-1]);
        let z = &un.slice(s![..-1, 1..]);
        let rhs = x - (c * dt / dx * (x - y)) - (c * dt / dy * (x - z));
        u.slice_mut(s![1.., 1..]).assign(&rhs);
        u.slice_mut(s![0, ..]).fill(1.);
        u.slice_mut(s![-1, ..]).fill(1.);
        u.slice_mut(s![.., 0]).fill(1.);
        u.slice_mut(s![.., -1]).fill(1.);
    }
    write_npy("array.npy", &u).unwrap();
}
