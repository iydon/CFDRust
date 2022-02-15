use ndarray::prelude::*;
use ndarray_npy::write_npy;

fn main() {
    // \pdv{u}{t} + u\pdv{u}{x} + v\pdv{u}{y} = ν(\pdv[2]{u}{x} + \pdv[2]{u}{y})
    // \pdv{v}{t} + u\pdv{v}{x} + v\pdv{v}{y} = ν(\pdv[2]{v}{x} + \pdv[2]{v}{y})
    let nx = 41;
    let ny = 41;
    let nt = 120;
    let nu = 0.01;
    let sigma = 0.0009;

    let dx = 2. / (nx as f64 - 1.);
    let dy = 2. / (ny as f64 - 1.);
    let dt = sigma * dx * dy / nu;
    let mut u = Array::ones((ny, nx));
    let mut v = Array::ones((ny, nx));
    let begin = (0.5 / dx) as usize;
    let end = (1. / dx + 1.) as usize;
    u.slice_mut(s![begin..end, begin..end]).fill(2.);  // nx == ny
    v.slice_mut(s![begin..end, begin..end]).fill(2.);  // nx == ny

    for _ in 0..nt + 1 {
        let un = u.clone();
        let vn = v.clone();
        let u0 = &un.slice(s![1..-1, 1..-1]);
        let u1 = &un.slice(s![1..-1, 2..]);
        let u2 = &un.slice(s![1..-1, ..-2]);
        let u3 = &un.slice(s![2.., 1..-1]);
        let u4 = &un.slice(s![..-2, 1..-1]);
        let v0 = &vn.slice(s![1..-1, 1..-1]);
        let v1 = &vn.slice(s![1..-1, 2..]);
        let v2 = &vn.slice(s![1..-1, ..-2]);
        let v3 = &vn.slice(s![2.., 1..-1]);
        let v4 = &vn.slice(s![..-2, 1..-1]);
        let rhs = u0 - dt / dx * u0 * (u0 - u2) - dt / dy * v0 * (u0 - u4)
            + nu * dt / dx.powi(2) * (u1 - 2. * u0 + u2)
            + nu * dt / dy.powi(2) * (u3 - 2. * u0 + u4);
        u.slice_mut(s![1..-1, 1..-1]).assign(&rhs);
        let rhs = v0 - dt / dx * u0 * (v0 - v2) - dt / dy * v0 * (v0 - v4)
            + nu * dt / dx.powi(2) * (v1 - 2. * v0 + v2)
            + nu * dt / dy.powi(2) * (v3 - 2. * v0 + v4);
        v.slice_mut(s![1..-1, 1..-1]).assign(&rhs);
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
