// https://github.com/barbagroup/CFDPython/blob/master/lessons/14_Step_11.ipynb
use ndarray::prelude::*;

use crate::macros::{assign, fill, i};

pub fn default() -> (
    Array2<f64>,
    Array2<f64>,
    Array2<f64>,
    usize,
    usize,
    f64,
    f64,
    f64,
    f64,
    f64,
) {
    let nx = 41;
    let ny = 41;
    let nt = 500;
    let nit = 50;
    let dt = 0.001;

    let rho = 1.;
    let nu = 0.1;

    let dx = 2. / (nx as f64 - 1.);
    let dy = 2. / (ny as f64 - 1.);

    let u = Array::zeros((ny, nx));
    let v = Array::zeros((ny, nx));
    let p = Array::zeros((ny, nx));

    return (u, v, p, nt, nit, dx, dy, dt, rho, nu);
}

pub fn solve(
    u: &mut Array2<f64>,
    v: &mut Array2<f64>,
    p: &mut Array2<f64>,
    nt: usize,
    nit: usize,
    dx: f64,
    dy: f64,
    dt: f64,
    rho: f64,
    nu: f64,
) {
    let mut b = Array::zeros(u.raw_dim());
    for _ in 0..nt {
        let un = u.clone();
        let vn = v.clone();
        set_b(&mut b, &u, &v, dx, dy, dt, rho);
        set_pressure_poisson(p, &b, dx, dy, nit);

        assign!(
            u[1..-1, 1..-1] =  i!(un[1..-1, 1..-1])
                - i!(un[1..-1, 1..-1]) * dt / dx * (i!(un[1..-1, 1..-1]) - i!(un[1..-1, ..-2]))
                - i!(vn[1..-1, 1..-1]) * dt / dy * (i!(un[1..-1, 1..-1]) - i!(un[..-2, 1..-1]))
                - dt / (2. * rho * dx) * (i!(p[1..-1, 2..]) - i!(p[1..-1, ..-2]))
                + nu * (dt / dx.powi(2)
                    * (i!(un[1..-1, 2..]) - 2. * i!(un[1..-1, 1..-1]) + i!(un[1..-1, ..-2]))
                    + dt / dy.powi(2)
                        * (i!(un[2.., 1..-1]) - 2. * i!(un[1..-1, 1..-1]) + i!(un[..-2, 1..-1])))
        );
        assign!(
            v[1..-1, 1..-1] = i!(vn[1..-1, 1..-1])
                - i!(un[1..-1, 1..-1]) * dt / dx * (i!(vn[1..-1, 1..-1]) - i!(vn[1..-1, ..-2]))
                - i!(vn[1..-1, 1..-1]) * dt / dy * (i!(vn[1..-1, 1..-1]) - i!(vn[..-2, 1..-1]))
                - dt / (2. * rho * dy) * (i!(p[2.., 1..-1]) - i!(p[..-2, 1..-1]))
                + nu * (dt / dx.powi(2)
                    * (i!(vn[1..-1, 2..]) - 2. * i!(vn[1..-1, 1..-1]) + i!(vn[1..-1, ..-2]))
                    + dt / dy.powi(2)
                        * (i!(vn[2.., 1..-1]) - 2. * i!(vn[1..-1, 1..-1]) + i!(vn[..-2, 1..-1])))
        );

        fill!(u[0, ..] = 0.);
        fill!(u[-1, ..] = 1.);
        fill!(u[.., 0] = 0.);
        fill!(u[.., -1] = 0.);
        fill!(v[0, ..] = 0.);
        fill!(v[-1, ..] = 0.);
        fill!(v[.., 0] = 0.);
        fill!(v[.., -1] = 0.);
    }
}

fn set_b(
    b: &mut Array2<f64>,
    u: &Array2<f64>,
    v: &Array2<f64>,
    dx: f64,
    dy: f64,
    dt: f64,
    rho: f64,
) {
    assign!(
        b[1..-1, 1..-1] = rho
            * (1. / dt
                * ((i!(u[1..-1, 2..]) - i!(u[1..-1, ..-2])) / (2. * dx)
                    + (i!(v[2.., 1..-1]) - i!(v[..-2, 1..-1])) / (2. * dy))
                - ((i!(u[1..-1, 2..]) - i!(u[1..-1, ..-2])) / (2. * dx)).mapv(|x| x.powi(2))
                - 2. * ((i!(u[2.., 1..-1]) - i!(u[..-2, 1..-1])) / (2. * dy)
                    * (i!(v[1..-1, 2..]) - i!(v[1..-1, ..-2]))
                    / (2. * dx))
                - ((i!(v[2.., 1..-1]) - i!(v[..-2, 1..-1])) / (2. * dy)).mapv(|x| x.powi(2)))
    );
}

fn set_pressure_poisson(p: &mut Array2<f64>, b: &Array2<f64>, dx: f64, dy: f64, nit: usize) {
    for _ in 0..nit {
        let pn = p.clone();

        assign!(
            p[1..-1, 1..-1] = ((i!(pn[1..-1, 2..]) + i!(pn[1..-1, ..-2])) * dy.powi(2)
                + (i!(pn[2.., 1..-1]) + i!(pn[..-2, 1..-1])) * dx.powi(2))
                / (2. * (dx.powi(2) + dy.powi(2)))
                - dx.powi(2) * dy.powi(2) / (2. * (dx.powi(2) + dy.powi(2))) * i!(b[1..-1, 1..-1])
        );
    }
}
