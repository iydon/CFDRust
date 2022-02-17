use ndarray::prelude::*;

use crate::macros::{assign, fill, i};

pub fn default() -> (
    Array2<f64>,
    Array2<f64>,
    Array2<f64>,
    usize,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
) {
    let nx = 41;
    let ny = 41;
    let nit = 50;
    let dt = 0.01;

    let rho = 1.;
    let nu = 0.1;
    let f = 1.;
    let eps = 0.001;

    let dx = 2. / (nx as f64 - 1.);
    let dy = 2. / (ny as f64 - 1.);

    let mut u = Array::zeros((ny, nx));
    let mut v = Array::zeros((ny, nx));
    let mut p = Array::ones((ny, nx));

    return (u, v, p, nit, dx, dy, dt, rho, nu, f, eps);
}

pub fn solve(
    u: &mut Array2<f64>,
    v: &mut Array2<f64>,
    p: &mut Array2<f64>,
    nit: usize,
    dx: f64,
    dy: f64,
    dt: f64,
    rho: f64,
    nu: f64,
    f: f64,
    eps: f64,
) {
    let mut udiff = 1.;
    let mut stepcount = 0;
    while udiff > 0.001 {
        let un = u.clone();
        let vn = v.clone();
        let b = make_b(&u, &v, dx, dy, dt, rho);
        set_pressure_poisson(p, b, nit, dx, dy);

        assign!(
            u[1..-1, 1..-1] = i!(un[1..-1, 1..-1])
                - i!(un[1..-1, 1..-1]) * dt / dx * (i!(un[1..-1, 1..-1]) - i!(un[1..-1, 0..-2]))
                - i!(vn[1..-1, 1..-1]) * dt / dy * (i!(un[1..-1, 1..-1]) - i!(un[0..-2, 1..-1]))
                - dt / (2. * rho * dx) * (i!(p[1..-1, 2..]) - i!(p[1..-1, 0..-2]))
                + nu * (dt / dx.powi(2)
                    * (i!(un[1..-1, 2..]) - 2. * i!(un[1..-1, 1..-1]) + i!(un[1..-1, 0..-2]))
                    + dt / dy.powi(2)
                        * (i!(un[2.., 1..-1]) - 2. * i!(un[1..-1, 1..-1]) + i!(un[0..-2, 1..-1])))
                + f * dt
        );
        assign!(
            v[1..-1, 1..-1] = i!(vn[1..-1, 1..-1])
                - i!(un[1..-1, 1..-1]) * dt / dx * (i!(vn[1..-1, 1..-1]) - i!(vn[1..-1, 0..-2]))
                - i!(vn[1..-1, 1..-1]) * dt / dy * (i!(vn[1..-1, 1..-1]) - i!(vn[0..-2, 1..-1]))
                - dt / (2. * rho * dy) * (i!(p[2.., 1..-1]) - i!(p[0..-2, 1..-1]))
                + nu * (dt / dx.powi(2)
                    * (i!(vn[1..-1, 2..]) - 2. * i!(vn[1..-1, 1..-1]) + i!(vn[1..-1, 0..-2]))
                    + dt / dy.powi(2)
                        * (i!(vn[2.., 1..-1]) - 2. * i!(vn[1..-1, 1..-1]) + i!(vn[0..-2, 1..-1])))
        );

        // Periodic BC u @ x = 2
        assign!(
            u[1..-1, -1] = i!(un[1..-1, -1])
                - i!(un[1..-1, -1]) * dt / dx * (i!(un[1..-1, -1]) - i!(un[1..-1, -2]))
                - i!(vn[1..-1, -1]) * dt / dy * (i!(un[1..-1, -1]) - i!(un[0..-2, -1]))
                - dt / (2. * rho * dx) * (i!(p[1..-1, 0]) - i!(p[1..-1, -2]))
                + nu * (dt / dx.powi(2)
                    * (i!(un[1..-1, 0]) - 2. * i!(un[1..-1,-1]) + i!(un[1..-1, -2]))
                    + dt / dy.powi(2) * (i!(un[2.., -1]) - 2. * i!(un[1..-1, -1]) + i!(un[0..-2, -1])))
                + f * dt
        );
        // Periodic BC u @ x = 0
        assign!(
            u[1..-1, 0] = i!(un[1..-1, 0])
                - i!(un[1..-1, 0]) * dt / dx * (i!(un[1..-1, 0]) - i!(un[1..-1, -1]))
                - i!(vn[1..-1, 0]) * dt / dy * (i!(un[1..-1, 0]) - i!(un[0..-2, 0]))
                - dt / (2. * rho * dx) * (i!(p[1..-1, 1]) - i!(p[1..-1, -1]))
                + nu * (dt / dx.powi(2)
                    * (i!(un[1..-1, 1]) - 2. * i!(un[1..-1, 0]) + i!(un[1..-1, -1]))
                    + dt / dy.powi(2) * (i!(un[2.., 0]) - 2. * i!(un[1..-1, 0]) + i!(un[0..-2, 0])))
                + f * dt
        );
        // Periodic BC v @ x = 2
        assign!(
            v[1..-1, -1] = i!(vn[1..-1, -1])
                - i!(un[1..-1, -1]) * dt / dx * (i!(vn[1..-1, -1]) - i!(vn[1..-1, -2]))
                - i!(vn[1..-1, -1]) * dt / dy * (i!(vn[1..-1, -1]) - i!(vn[0..-2, -1]))
                - dt / (2. * rho * dy) * (i!(p[2.., -1]) - i!(p[0..-2, -1]))
                + nu * (dt / dx.powi(2)
                    * (i!(vn[1..-1, 0]) - 2. * i!(vn[1..-1, -1]) + i!(vn[1..-1, -2]))
                    + dt / dy.powi(2) * (i!(vn[2.., -1]) - 2. * i!(vn[1..-1, -1]) + i!(vn[0..-2, -1])))
        );
        // Periodic BC v @ x = 0
        assign!(
            v[1..-1, 0] = i!(vn[1..-1, 0])
                - i!(un[1..-1, 0]) * dt / dx * (i!(vn[1..-1, 0]) - i!(vn[1..-1, -1]))
                - i!(vn[1..-1, 0]) * dt / dy * (i!(vn[1..-1, 0]) - i!(vn[0..-2, 0]))
                - dt / (2. * rho * dy) * (i!(p[2.., 0]) - i!(p[0..-2, 0]))
                + nu * (dt / dx.powi(2)
                    * (i!(vn[1..-1, 1]) - 2. * i!(vn[1..-1, 0]) + i!(vn[1..-1, -1]))
                    + dt / dy.powi(2) * (i!(vn[2.., 0]) - 2. * i!(vn[1..-1, 0]) + i!(vn[0..-2, 0])))
        );
        // Wall BC: u,v = 0 @ y = 0, 2
        fill!(u[0, ..] = 0.);
        fill!(u[-1, ..] = 0.);
        fill!(v[0, ..] = 0.);
        fill!(v[-1, ..] = 0.);

        udiff = (u.sum() - un.sum()).abs() / u.sum();
        stepcount += 1;
    }
    println!("{}", stepcount);
}

fn make_b(u: &Array2<f64>, v: &Array2<f64>, dx: f64, dy: f64, dt: f64, rho: f64) -> Array2<f64> {
    let mut b = Array::zeros(u.raw_dim());
    assign!(
        b[1..-1, 1..-1] = rho
            * (1. / dt
                * ((i!(u[1..-1, 2..]) - i!(u[1..-1, 0..-2])) / (2. * dx)
                    + (i!(v[2.., 1..-1]) - i!(v[0..-2, 1..-1])) / (2. * dy))
                - ((i!(u[1..-1, 2..]) - i!(u[1..-1, 0..-2])) / (2. * dx)).mapv(|x| x.powi(2))
                - 2. * ((i!(u[2.., 1..-1]) - i!(u[0..-2, 1..-1])) / (2. * dy)
                    * (i!(v[1..-1, 2..]) - i!(v[1..-1, 0..-2]))
                    / (2. * dx))
                - ((i!(v[2.., 1..-1]) - i!(v[0..-2, 1..-1])) / (2. * dy)).mapv(|x| x.powi(2)))
    );
    // Periodic BC Pressure @ x = 2
    assign!(
        b[1..-1, -1] = rho
            * (1. / dt
                * ((i!(u[1..-1, 0]) - i!(u[1..-1,-2])) / (2. * dx)
                    + (i!(v[2.., -1]) - i!(v[0..-2, -1])) / (2. * dy))
                - ((i!(u[1..-1, 0]) - i!(u[1..-1, -2])) / (2. * dx)).mapv(|x| x.powi(2))
                - 2. * ((i!(u[2.., -1]) - i!(u[0..-2, -1])) / (2. * dy)
                    * (i!(v[1..-1, 0]) - i!(v[1..-1, -2]))
                    / (2. * dx))
                - ((i!(v[2.., -1]) - i!(v[0..-2, -1])) / (2. * dy)).mapv(|x| x.powi(2)))
    );
    // Periodic BC Pressure @ x = 0
    assign!(
        b[1..-1, 0] = rho
            * (1. / dt
                * ((i!(u[1..-1, 1]) - i!(u[1..-1, -1])) / (2. * dx)
                    + (i!(v[2.., 0]) - i!(v[0..-2, 0])) / (2. * dy))
                - ((i!(u[1..-1, 1]) - i!(u[1..-1, -1])) / (2. * dx)).mapv(|x| x.powi(2))
                - 2. * ((i!(u[2.., 0]) - i!(u[0..-2, 0])) / (2. * dy)
                    * (i!(v[1..-1, 1]) - i!(v[1..-1, -1]))
                    / (2. * dx))
                - ((i!(v[2.., 0]) - i!(v[0..-2, 0])) / (2. * dy)).mapv(|x| x.powi(2)))
    );
    return b;
}

fn set_pressure_poisson(p: &mut Array2<f64>, b: Array2<f64>, nit: usize, dx: f64, dy: f64) {
    for _ in 0..nit {
        let pn = p.clone();

        assign!(
            p[1..-1, 1..-1] = ((i!(pn[1..-1, 2..]) + i!(pn[1..-1, 0..-2])) * dy.powi(2)
                + (i!(pn[2.., 1..-1]) + i!(pn[0..-2, 1..-1])) * dx.powi(2))
                / (2. * (dx.powi(2) + dy.powi(2)))
                - dx.powi(2) * dy.powi(2) / (2. * (dx.powi(2) + dy.powi(2))) * i!(b[1..-1, 1..-1])
        );
        // Periodic BC Pressure @ x = 2
        assign!(
            p[1..-1, -1] = ((i!(pn[1..-1, 0]) + i!(pn[1..-1, -2])) * dy.powi(2)
                + (i!(pn[2.., -1]) + i!(pn[0..-2, -1])) * dx.powi(2))
                / (2. * (dx.powi(2) + dy.powi(2)))
                - dx.powi(2) * dy.powi(2) / (2. * (dx.powi(2) + dy.powi(2))) * i!(b[1..-1, -1])
        );
        // Periodic BC Pressure @ x = 0
        assign!(
            p[1..-1, 0] = ((i!(pn[1..-1, 1]) + i!(pn[1..-1, -1])) * dy.powi(2)
                + (i!(pn[2.., 0]) + i!(pn[0..-2, 0])) * dx.powi(2))
                / (2. * (dx.powi(2) + dy.powi(2)))
                - dx.powi(2) * dy.powi(2) / (2. * (dx.powi(2) + dy.powi(2))) * i!(b[1..-1, 0])
        );
        // Wall boundary conditions, pressure
        let pn = i!(p[1, ..]).to_owned();
        assign!(p[0, ..] = pn); // dp/dy = 0 @ y = 0
        let pn = i!(p[-2, ..]).to_owned();
        assign!(p[-1, ..] = pn); // dp/dy = 0 @ y = 2
    }
}
