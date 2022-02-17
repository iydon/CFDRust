use ndarray::prelude::*;

pub fn default() -> (Array1<f64>, usize, f64, f64) {
    let nx = 41;
    let nt = 25;
    let dt = 0.025;

    let dx = 2. / (nx as f64 - 1.);

    let mut u = Array::ones(nx);
    let begin = (0.5 / dx) as usize;
    let end = (1. / dx + 1.) as usize;
    u.slice_mut(s![begin..end]).fill(2.);

    return (u, nt, dx, dt);
}

pub fn solve(u: &mut Array1<f64>, nt: usize, dx: f64, dt: f64) {
    let nx = u.shape()[0];
    for _ in 0..nt {
        let un = u.clone();
        for ith in 1..nx {
            u[ith] = un[ith] - un[ith] * dt / dx * (un[ith] - un[ith - 1]);
        }
    }
}
