use ndarray::prelude::*;

fn main() {
    // \pdv{u}{t} + u\pdv{u}{x} = 0
    let nx = 85;
    let nt = 25;
    let sigma = 0.5;
    let c = 1.0;

    let dx = 2. / (nx as f64 - 1.);
    let dt = sigma * dx;
    let mut u = Array::ones(nx);
    let begin = (0.5 / dx) as usize;
    let end = (1. / dx + 1.) as usize;
    u.slice_mut(s![begin..end]).fill(2.);

    for _ in 0..nt {
        let un = u.clone();
        for ith in 1..nx {
            u[ith] = un[ith] - c * dt / dx * (un[ith] - un[ith - 1]);
        }
    }
    println!("{}", u);
}