use ndarray::prelude::*;
use ndarray_npy::write_npy;

macro_rules! i {
    ($v:ident[$($x:expr),*]) => {
        &$v.slice(s![$($x),*])
    };
}

fn build_up_b(
    u: &Array2<f64>,
    v: &Array2<f64>,
    rho: f64,
    dt: f64,
    dx: f64,
    dy: f64,
) -> Array2<f64> {
    let mut b = Array::zeros(u.raw_dim());
    let rhs = rho
        * (1. / dt
            * ((i!(u[1..-1, 2..]) - i!(u[1..-1, 0..-2])) / (2. * dx)
                + (i!(v[2.., 1..-1]) - i!(v[0..-2, 1..-1])) / (2. * dy))
            - ((i!(u[1..-1, 2..]) - i!(u[1..-1, 0..-2])) / (2. * dx)).mapv(|x| x.powi(2))
            - 2. * ((i!(u[2.., 1..-1]) - i!(u[0..-2, 1..-1])) / (2. * dy)
                * (i!(v[1..-1, 2..]) - i!(v[1..-1, 0..-2]))
                / (2. * dx))
            - ((i!(v[2.., 1..-1]) - i!(v[0..-2, 1..-1])) / (2. * dy)).mapv(|x| x.powi(2)));
    b.slice_mut(s![1..-1, 1..-1]).assign(&rhs);
    // Periodic BC Pressure @ x = 2
    let rhs = rho
        * (1. / dt
            * ((i!(u[1..-1, 0]) - i!(u[1..-1,-2])) / (2. * dx)
                + (i!(v[2.., -1]) - i!(v[0..-2, -1])) / (2. * dy))
            - ((i!(u[1..-1, 0]) - i!(u[1..-1, -2])) / (2. * dx)).mapv(|x| x.powi(2))
            - 2. * ((i!(u[2.., -1]) - i!(u[0..-2, -1])) / (2. * dy)
                * (i!(v[1..-1, 0]) - i!(v[1..-1, -2]))
                / (2. * dx))
            - ((i!(v[2.., -1]) - i!(v[0..-2, -1])) / (2. * dy)).mapv(|x| x.powi(2)));
    b.slice_mut(s![1..-1, -1]).assign(&rhs);
    // Periodic BC Pressure @ x = 0
    let rhs = rho
        * (1. / dt
            * ((i!(u[1..-1, 1]) - i!(u[1..-1, -1])) / (2. * dx)
                + (i!(v[2.., 0]) - i!(v[0..-2, 0])) / (2. * dy))
            - ((i!(u[1..-1, 1]) - i!(u[1..-1, -1])) / (2. * dx)).mapv(|x| x.powi(2))
            - 2. * ((i!(u[2.., 0]) - i!(u[0..-2, 0])) / (2. * dy)
                * (i!(v[1..-1, 1]) - i!(v[1..-1, -1]))
                / (2. * dx))
            - ((i!(v[2.., 0]) - i!(v[0..-2, 0])) / (2. * dy)).mapv(|x| x.powi(2)));
    b.slice_mut(s![1..-1, 0]).assign(&rhs);
    return b;
}

fn pressure_poisson(p: &mut Array2<f64>, b: Array2<f64>, dx: f64, dy: f64, nit: usize) {
    for _ in 0..nit {
        let pn = p.clone();
        let rhs = ((i!(pn[1..-1, 2..]) + i!(pn[1..-1, 0..-2])) * dy.powi(2)
            + (i!(pn[2.., 1..-1]) + i!(pn[0..-2, 1..-1])) * dx.powi(2))
            / (2. * (dx.powi(2) + dy.powi(2)))
            - dx.powi(2) * dy.powi(2) / (2. * (dx.powi(2) + dy.powi(2))) * i!(b[1..-1, 1..-1]);
        p.slice_mut(s![1..-1, 1..-1]).assign(&rhs);
        // Periodic BC Pressure @ x = 2
        let rhs = ((i!(pn[1..-1, 0]) + i!(pn[1..-1, -2])) * dy.powi(2)
            + (i!(pn[2.., -1]) + i!(pn[0..-2, -1])) * dx.powi(2))
            / (2. * (dx.powi(2) + dy.powi(2)))
            - dx.powi(2) * dy.powi(2) / (2. * (dx.powi(2) + dy.powi(2))) * i!(b[1..-1, -1]);
        p.slice_mut(s![1..-1, -1]).assign(&rhs);
        // Periodic BC Pressure @ x = 0
        let rhs = ((i!(pn[1..-1, 1]) + i!(pn[1..-1, -1])) * dy.powi(2)
            + (i!(pn[2.., 0]) + i!(pn[0..-2, 0])) * dx.powi(2))
            / (2. * (dx.powi(2) + dy.powi(2)))
            - dx.powi(2) * dy.powi(2) / (2. * (dx.powi(2) + dy.powi(2))) * i!(b[1..-1, 0]);
        p.slice_mut(s![1..-1, 0]).assign(&rhs);
        // Wall boundary conditions, pressure
        let p1 = p.slice(s![1, ..]).to_owned();
        p.slice_mut(s![0, ..]).assign(&p1); // dp/dy = 0 @ y = 0
        let p2 = p.slice(s![-2, ..]).to_owned();
        p.slice_mut(s![-1, ..]).assign(&p2); // dp/dy = 0 @ y = 2
    }
}

fn main() {
    let nx = 41;
    let ny = 41;
    let nit = 50;

    let rho = 1.;
    let nu = 0.1;
    let f = 1.;
    let dt = 0.01;

    let dx = 2. / (nx as f64 - 1.);
    let dy = 2. / (ny as f64 - 1.);
    let mut u = Array::zeros((ny, nx));
    let mut v = Array::zeros((ny, nx));
    let mut p = Array::zeros((ny, nx));

    let mut udiff = 1.;
    let mut stepcount = 0;
    while udiff > 0.001 {
        let un = u.clone();
        let vn = v.clone();
        let b = build_up_b(&u, &v, rho, dt, dx, dy);
        pressure_poisson(&mut p, b, dx, dy, nit);

        let rhs = i!(un[1..-1, 1..-1])
            - i!(un[1..-1, 1..-1]) * dt / dx * (i!(un[1..-1, 1..-1]) - i!(un[1..-1, 0..-2]))
            - i!(vn[1..-1, 1..-1]) * dt / dy * (i!(un[1..-1, 1..-1]) - i!(un[0..-2, 1..-1]))
            - dt / (2. * rho * dx) * (i!(p[1..-1, 2..]) - i!(p[1..-1, 0..-2]))
            + nu * (dt / dx.powi(2)
                * (i!(un[1..-1, 2..]) - 2. * i!(un[1..-1, 1..-1]) + i!(un[1..-1, 0..-2]))
                + dt / dy.powi(2)
                    * (i!(un[2.., 1..-1]) - 2. * i!(un[1..-1, 1..-1]) + i!(un[0..-2, 1..-1])))
            + f * dt;
        u.slice_mut(s![1..-1, 1..-1]).assign(&rhs);
        let rhs = i!(vn[1..-1, 1..-1])
            - i!(un[1..-1, 1..-1]) * dt / dx * (i!(vn[1..-1, 1..-1]) - i!(vn[1..-1, 0..-2]))
            - i!(vn[1..-1, 1..-1]) * dt / dy * (i!(vn[1..-1, 1..-1]) - i!(vn[0..-2, 1..-1]))
            - dt / (2. * rho * dy) * (i!(p[2.., 1..-1]) - i!(p[0..-2, 1..-1]))
            + nu * (dt / dx.powi(2)
                * (i!(vn[1..-1, 2..]) - 2. * i!(vn[1..-1, 1..-1]) + i!(vn[1..-1, 0..-2]))
                + dt / dy.powi(2)
                    * (i!(vn[2.., 1..-1]) - 2. * i!(vn[1..-1, 1..-1]) + i!(vn[0..-2, 1..-1])));
        v.slice_mut(s![1..-1, 1..-1]).assign(&rhs);
        // Periodic BC u @ x = 2
        let rhs = i!(un[1..-1, -1])
            - i!(un[1..-1, -1]) * dt / dx * (i!(un[1..-1, -1]) - i!(un[1..-1, -2]))
            - i!(vn[1..-1, -1]) * dt / dy * (i!(un[1..-1, -1]) - i!(un[0..-2, -1]))
            - dt / (2. * rho * dx) * (i!(p[1..-1, 0]) - i!(p[1..-1, -2]))
            + nu * (dt / dx.powi(2)
                * (i!(un[1..-1, 0]) - 2. * i!(un[1..-1,-1]) + i!(un[1..-1, -2]))
                + dt / dy.powi(2) * (i!(un[2.., -1]) - 2. * i!(un[1..-1, -1]) + i!(un[0..-2, -1])))
            + f * dt;
        u.slice_mut(s![1..-1, -1]).assign(&rhs);
        // Periodic BC u @ x = 0
        let rhs = i!(un[1..-1, 0])
            - i!(un[1..-1, 0]) * dt / dx * (i!(un[1..-1, 0]) - i!(un[1..-1, -1]))
            - i!(vn[1..-1, 0]) * dt / dy * (i!(un[1..-1, 0]) - i!(un[0..-2, 0]))
            - dt / (2. * rho * dx) * (i!(p[1..-1, 1]) - i!(p[1..-1, -1]))
            + nu * (dt / dx.powi(2)
                * (i!(un[1..-1, 1]) - 2. * i!(un[1..-1, 0]) + i!(un[1..-1, -1]))
                + dt / dy.powi(2) * (i!(un[2.., 0]) - 2. * i!(un[1..-1, 0]) + i!(un[0..-2, 0])))
            + f * dt;
        u.slice_mut(s![1..-1, 0]).assign(&rhs);
        // Periodic BC v @ x = 2
        let rhs = i!(vn[1..-1, -1])
            - i!(un[1..-1, -1]) * dt / dx * (i!(vn[1..-1, -1]) - i!(vn[1..-1, -2]))
            - i!(vn[1..-1, -1]) * dt / dy * (i!(vn[1..-1, -1]) - i!(vn[0..-2, -1]))
            - dt / (2. * rho * dy) * (i!(p[2.., -1]) - i!(p[0..-2, -1]))
            + nu * (dt / dx.powi(2)
                * (i!(vn[1..-1, 0]) - 2. * i!(vn[1..-1, -1]) + i!(vn[1..-1, -2]))
                + dt / dy.powi(2) * (i!(vn[2.., -1]) - 2. * i!(vn[1..-1, -1]) + i!(vn[0..-2, -1])));
        v.slice_mut(s![1..-1, -1]).assign(&rhs);
        // Periodic BC v @ x = 0
        let rhs = i!(vn[1..-1, 0])
            - i!(un[1..-1, 0]) * dt / dx * (i!(vn[1..-1, 0]) - i!(vn[1..-1, -1]))
            - i!(vn[1..-1, 0]) * dt / dy * (i!(vn[1..-1, 0]) - i!(vn[0..-2, 0]))
            - dt / (2. * rho * dy) * (i!(p[2.., 0]) - i!(p[0..-2, 0]))
            + nu * (dt / dx.powi(2)
                * (i!(vn[1..-1, 1]) - 2. * i!(vn[1..-1, 0]) + i!(vn[1..-1, -1]))
                + dt / dy.powi(2) * (i!(vn[2.., 0]) - 2. * i!(vn[1..-1, 0]) + i!(vn[0..-2, 0])));
        v.slice_mut(s![1..-1, 0]).assign(&rhs);
        // Wall BC: u,v = 0 @ y = 0, 2
        u.slice_mut(s![0, ..]).fill(0.);
        u.slice_mut(s![-1, ..]).fill(0.);
        v.slice_mut(s![0, ..]).fill(0.);
        v.slice_mut(s![0, ..]).fill(0.);

        udiff = (u.sum() - un.sum()).abs() / u.sum();
        stepcount += 1;
    }
    println!("{}", stepcount);

    write_npy("array_u.npy", &u).unwrap();
    write_npy("array_v.npy", &v).unwrap();
    write_npy("array_p.npy", &p).unwrap();
}
