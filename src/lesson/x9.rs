// https://github.com/barbagroup/CFDPython/blob/master/lessons/09_Step_7.ipynb
use ndarray::prelude::*;

use crate::macros::{assign, fill, i};

pub fn default() -> (Array2<f64>, usize, f64, f64, f64, f64) {
    let nx = 31;
    let ny = 31;
    let nt = 17;

    let nu = 0.05;
    let sigma = 0.25;

    let dx = 2. / (nx as f64 - 1.);
    let dy = 2. / (ny as f64 - 1.);
    let dt = sigma * dx * dy / nu;

    let mut u = Array::ones((ny, nx));
    let x_begin = (0.5 / dx) as usize;
    let x_end = (1. / dx + 1.) as usize;
    let y_begin = (0.5 / dy) as usize;
    let y_end = (1. / dy + 1.) as usize;
    fill!(u[x_begin..x_end, y_begin..y_end] = 2.);

    return (u, nt, dx, dy, dt, nu);
}

pub fn solve(u: &mut Array2<f64>, nt: usize, dx: f64, dy: f64, dt: f64, nu: f64) {
    for _ in 0..nt + 1 {
        let un = u.clone();

        assign!(
            u[1..-1, 1..-1] = i!(un[1..-1, 1..-1])
                + nu * dt / dx.powi(2)
                    * (i!(un[1..-1, 2..]) - 2. * i!(un[1..-1, 1..-1]) + i!(un[1..-1, ..-2]))
                + nu * dt / dy.powi(2)
                    * (i!(un[2.., 1..-1]) - 2. * i!(un[1..-1, 1..-1]) + i!(un[..-2, 1..-1]))
        );

        fill!(u[0, ..] = 1.);
        fill!(u[-1, ..] = 1.);
        fill!(u[.., 0] = 1.);
        fill!(u[.., -1] = 1.);
    }
}
