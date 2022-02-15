use ndarray::prelude::*;
use ndarray_npy::write_npy;

fn build_up_b(
    b: &mut Array2<f64>,
    u: &Array2<f64>,
    v: &Array2<f64>,
    rho: f64,
    dt: f64,
    dx: f64,
    dy: f64,
) {
    let u1 = &u.slice(s![1..-1, 2..]);
    let u2 = &u.slice(s![1..-1, ..-2]);
    let u3 = &u.slice(s![2.., 1..-1]);
    let u4 = &u.slice(s![..-2, 1..-1]);
    let v1 = &v.slice(s![2.., 1..-1]);
    let v2 = &v.slice(s![..-2, 1..-1]);
    let v3 = &v.slice(s![1..-1, 2..]);
    let v4 = &v.slice(s![1..-1, ..-2]);
    let rhs = rho
        * (1. / dt * ((u1 - u2) / (2. * dx) + (v1 - v2) / (2. * dy))
            - ((u1 - u2) / (2. * dx)).mapv(|x| x.powi(2))
            - 2. * ((u3 - u4) / (2. * dy) * (v3 - v4) / (2. * dx))
            - ((v1 - v2) / (2. * dy)).mapv(|x| x.powi(2)));
    b.slice_mut(s![1..-1, 1..-1]).assign(&rhs);
}

fn pressure_poisson(p: &mut Array2<f64>, b: &Array2<f64>, dx: f64, dy: f64, nit: usize) {
    for _ in 0..nit {
        let pn = p.clone();
        let p1 = &pn.slice(s![1..-1, 2..]);
        let p2 = &pn.slice(s![1..-1, ..-2]);
        let p3 = &pn.slice(s![2.., 1..-1]);
        let p4 = &pn.slice(s![..-2, 1..-1]);
        let b0 = &b.slice(s![1..-1, 1..-1]);
        let rhs = ((p1 + p2) * dy.powi(2) + (p3 + p4) * dx.powi(2))
            / (2. * (dx.powi(2) + dy.powi(2)))
            - dx.powi(2) * dy.powi(2) / (2. * (dx.powi(2) + dy.powi(2))) * b0;
        p.slice_mut(s![1..-1, 1..-1]).assign(&rhs);
    }
}

fn cavity_flow(
    u: &mut Array2<f64>,
    v: &mut Array2<f64>,
    p: &mut Array2<f64>,
    rho: f64,
    nu: f64,
    dt: f64,
    nt: usize,
    nx: usize,
    ny: usize,
    nit: usize,
) {
    let dx = 2. / (nx as f64 - 1.);
    let dy = 2. / (ny as f64 - 1.);
    let mut b = Array::zeros((ny, nx));
    for _ in 0..nt {
        let un = u.clone();
        let vn = v.clone();
        build_up_b(&mut b, &u, &v, rho, dt, dx, dy);
        pressure_poisson(p, &b, dx, dy, nit);

        let u0 = &un.slice(s![1..-1, 1..-1]);
        let u1 = &un.slice(s![1..-1, ..-2]);
        let u2 = &un.slice(s![..-2, 1..-1]);
        let u3 = &un.slice(s![1..-1, 2..]);
        let u4 = &un.slice(s![2.., 1..-1]);
        let v0 = &vn.slice(s![1..-1, 1..-1]);
        let v1 = &vn.slice(s![1..-1, ..-2]);
        let v2 = &vn.slice(s![..-2, 1..-1]);
        let v3 = &vn.slice(s![1..-1, 2..]);
        let v4 = &vn.slice(s![2.., 1..-1]);
        let p1 = &p.slice(s![1..-1, 2..]);
        let p2 = &p.slice(s![1..-1, ..-2]);
        let p3 = &p.slice(s![2.., 1..-1]);
        let p4 = &p.slice(s![..-2, 1..-1]);
        let rhs = u0
            - u0 * dt / dx * (u0 - u1)
            - v0 * dt / dy * (u0 - u2)
            - dt / (2. * rho * dx) * (p1 - p2)
            + nu * (dt / dx.powi(2) * (u3 - 2. * u0 + u1) + dt / dy.powi(2) * (u4 - 2. * u0 + u2));
        u.slice_mut(s![1..-1, 1..-1]).assign(&rhs);
        let rhs = v0
            - u0 * dt / dx * (v0 - v1)
            - v0 * dt / dy * (v0 - v2)
            - dt / (2. * rho * dy) * (p3 - p4)
            + nu * (dt / dx.powi(2) * (v3 - 2. * v0 + v1) + dt / dy.powi(2) * (v4 - 2. * v0 + v2));
        v.slice_mut(s![1..-1, 1..-1]).assign(&rhs);

        u.slice_mut(s![0, ..]).fill(0.);
        u.slice_mut(s![.., 0]).fill(0.);
        u.slice_mut(s![.., -1]).fill(0.);
        u.slice_mut(s![-1, ..]).fill(1.);
        v.slice_mut(s![0, ..]).fill(0.);
        v.slice_mut(s![-1, ..]).fill(0.);
        v.slice_mut(s![.., 0]).fill(0.);
        v.slice_mut(s![.., -1]).fill(0.);
    }
}

fn main() {
    let nx = 41;
    let ny = 41;
    let nt = 500;
    let nit = 50;
    let rho = 1.;
    let nu = 0.1;
    let dt = 0.001;

    let mut u = Array::zeros((ny, nx));
    let mut v = Array::zeros((ny, nx));
    let mut p = Array::zeros((ny, nx));
    cavity_flow(&mut u, &mut v, &mut p, rho, nu, dt, nt, nx, ny, nit);
    write_npy("array_u.npy", &u).unwrap();
    write_npy("array_v.npy", &v).unwrap();
    write_npy("array_p.npy", &p).unwrap();
}
