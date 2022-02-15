use ndarray::prelude::*;
use ndarray_npy::write_npy;

fn main() {
    // \pdv{u}{t} + u\pdv{u}{x} + v\pdv{u}{y} = 0
    // \pdv{v}{t} + u\pdv{v}{x} + v\pdv{v}{y} = 0
    let nx = 101;
    let ny = 101;
    let nt = 80;
    let sigma = 0.2;

    let dx = 2. / (nx as f64 - 1.);
    let dy = 2. / (ny as f64 - 1.);
    let dt = sigma * dx;
    let mut u = Array::ones((nx, ny));
    let mut v = Array::ones((nx, ny));
    let begin = (0.5 / dx) as usize;
    let end = (1. / dx + 1.) as usize;
    u.slice_mut(s![begin..end, begin..end]).fill(2.);
    v.slice_mut(s![begin..end, begin..end]).fill(2.);

    for _ in 0..nt + 1 {
        let un = u.clone();
        let vn = v.clone();
        let u0 = &un.slice(s![1.., 1..]);
        let u1 = &un.slice(s![1.., ..-1]);
        let u2 = &un.slice(s![..-1, 1..]);
        let v0 = &vn.slice(s![1.., 1..]);
        let v1 = &vn.slice(s![1.., ..-1]);
        let v2 = &vn.slice(s![..-1, 1..]);
        let rhs = u0 - (u0 * dt / dx * (u0 - u1)) - (v0 * dt / dy * (u0 - u2));
        u.slice_mut(s![1.., 1..]).assign(&rhs);
        let rhs =v0 - (u0 * dt / dx * (v0 - v1)) - (v0 * dt / dy * (v0 - v2));
        v.slice_mut(s![1.., 1..]).assign(&rhs);
        u.slice_mut(s![0, ..]).fill(1.);
        u.slice_mut(s![-1, ..]).fill(1.);
        u.slice_mut(s![.., 0]).fill(1.);
        u.slice_mut(s![.., -1]).fill(1.);
        v.slice_mut(s![0, ..]).fill(1.);
        v.slice_mut(s![-1, ..]).fill(1.);
        v.slice_mut(s![.., 0]).fill(1.);
        v.slice_mut(s![.., -1]).fill(1.);
    }
    write_npy("array_u.npy", &u).unwrap();
    write_npy("array_v.npy", &u).unwrap();
}
