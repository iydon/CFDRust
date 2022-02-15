use ndarray::prelude::*;

fn main() {
    // \pdv{u}{t} = ν\pdv[2]{u}{x}
    let nx = 41;
    let nt = 20;
    let nu = 0.3;
    let sigma = 0.2;

    let dx = 2. / (nx as f64 - 1.);
    let dt = sigma * dx.powi(2) / nu;
    let mut u = Array::ones(nx);
    let begin = (0.5 / dx) as usize;
    let end = (1. / dx + 1.) as usize;
    u.slice_mut(s![begin..end]).fill(2.);

    for _ in 0..nt {
        let un = u.clone();
        for ith in 1..nx - 1 {
            u[ith] = un[ith] + nu * dt / dx.powi(2) * (un[ith + 1] - 2. * un[ith] + un[ith - 1]);
        }
    }
    println!("{}", u);
}
