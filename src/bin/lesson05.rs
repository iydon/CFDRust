use ndarray::prelude::*;
use std::f64::consts::PI;

fn ufunc(t: f64, x: f64, nu: f64) -> f64 {
    return -2.
        * nu
        * (-(-8. * t + 2. * x) * (-(-4. * t + x).powi(2) / (4. * nu * (t + 1.))).exp()
            / (4. * nu * (t + 1.))
            - (-8. * t + 2. * x - 4. * PI)
                * (-(-4. * t + x - 2. * PI).powi(2) / (4. * nu * (t + 1.))).exp()
                / (4. * nu * (t + 1.)))
        / ((-(-4. * t + x - 2. * PI).powi(2) / (4. * nu * (t + 1.))).exp()
            + (-(-4. * t + x).powi(2) / (4. * nu * (t + 1.))).exp())
        + 4.;
}

fn main() {
    // \pdv{u}{t} + u\pdv{u}{x} = ν\pdv[2]{u}{x}
    let nx = 101;
    let nt = 100;
    let nu = 0.07;

    let dx = 2. * PI / (nx as f64 - 1.);
    let dt = dx * nu;
    let mut u = Array::zeros(nx);
    for (ith, x0) in Array::linspace(0., 2. * PI, nx).into_iter().enumerate() {
        u[ith] = ufunc(0., x0, nu);
    }

    for _ in 0..nt {
        let un = u.clone();
        for ith in 1..nx - 1 {
            u[ith] = un[ith] - un[ith] * dt / dx * (un[ith] - un[ith - 1])
                + nu * dt / dx.powi(2) * (un[ith + 1] - 2. * un[ith] + un[ith - 1]);
        }
        u[0] = un[0] - un[0] * dt / dx * (un[0] - un[nx - 2])
            + nu * dt / dx.powi(2) * (un[1] - 2. * un[0] + un[nx - 2]);
        u[nx - 1] = u[0];
    }
    println!("{}", u);
}
