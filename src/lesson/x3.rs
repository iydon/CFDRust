// https://github.com/barbagroup/CFDPython/blob/master/lessons/03_CFL_Condition.ipynb
use ndarray::prelude::*;

use crate::macros::fill;

pub fn default() -> (Array1<f64>, usize, f64, f64, f64) {
    let nx = 85;
    let nt = 25;

    let c = 1.;
    let sigma = 0.5;

    let dx = 2. / (nx as f64 - 1.);
    let dt = sigma * dx;

    let mut u = Array::ones(nx);
    let begin = (0.5 / dx) as usize;
    let end = (1. / dx + 1.) as usize;
    fill!(u[begin..end] = 2.);

    return (u, nt, dx, dt, c);
}

pub fn solve(u: &mut Array1<f64>, nt: usize, dx: f64, dt: f64, c: f64) {
    let nx = u.shape()[0];
    for _ in 0..nt {
        let un = u.clone();

        for ith in 1..nx {
            u[ith] = un[ith] - c * dt / dx * (un[ith] - un[ith - 1]);
        }
    }
}
