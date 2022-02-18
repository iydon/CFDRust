// https://github.com/barbagroup/CFDPython/blob/master/lessons/04_Step_3.ipynb
use ndarray::prelude::*;

use crate::macros::fill;

pub fn default() -> (Array1<f64>, usize, f64, f64, f64) {
    let nx = 41;
    let nt = 20;

    let nu = 0.3;
    let sigma = 0.2;

    let dx = 2. / (nx as f64 - 1.);
    let dt = sigma * dx.powi(2) / nu;

    let mut u = Array::ones(nx);
    let begin = (0.5 / dx) as usize;
    let end = (1. / dx + 1.) as usize;
    fill!(u[begin..end] = 2.);

    return (u, nt, dx, dt, nu);
}

pub fn solve(u: &mut Array1<f64>, nt: usize, dx: f64, dt: f64, nu: f64) {
    let nx = u.shape()[0];
    for _ in 0..nt {
        let un = u.clone();

        for ith in 1..nx - 1 {
            u[ith] = un[ith] + nu * dt / dx.powi(2) * (un[ith + 1] - 2. * un[ith] + un[ith - 1]);
        }
    }
}
