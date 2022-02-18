// https://github.com/barbagroup/CFDPython/blob/master/lessons/08_Step_6.ipynb
use ndarray::prelude::*;

use crate::macros::{assign, fill, i};

pub fn default() -> (Array2<f64>, Array2<f64>, usize, f64, f64, f64) {
    let nx = 101;
    let ny = 101;
    let nt = 80;

    let sigma = 0.2;

    let dx = 2. / (nx as f64 - 1.);
    let dy = 2. / (ny as f64 - 1.);
    let dt = sigma * dx;

    let mut u = Array::ones((ny, nx));
    let mut v = Array::ones((ny, nx));
    let x_begin = (0.5 / dx) as usize;
    let x_end = (1. / dx + 1.) as usize;
    let y_begin = (0.5 / dy) as usize;
    let y_end = (1. / dy + 1.) as usize;
    fill!(u[x_begin..x_end, y_begin..y_end] = 2.);
    fill!(v[x_begin..x_end, y_begin..y_end] = 2.);

    return (u, v, nt, dx, dy, dt);
}

pub fn solve(u: &mut Array2<f64>, v: &mut Array2<f64>, nt: usize, dx: f64, dy: f64, dt: f64) {
    for _ in 0..nt + 1 {
        let un = u.clone();
        let vn = v.clone();

        assign!(
            u[1.., 1..] = i!(un[1.., 1..])
                - (i!(un[1.., 1..]) * dt / dx * (i!(un[1.., 1..]) - i!(un[1.., ..-1])))
                - (i!(vn[1.., 1..]) * dt / dy * (i!(un[1.., 1..]) - i!(un[..-1, 1..])))
        );
        assign!(
            v[1.., 1..] = i!(vn[1.., 1..])
                - (i!(un[1.., 1..]) * dt / dx * (i!(vn[1.., 1..]) - i!(vn[1.., ..-1])))
                - (i!(vn[1.., 1..]) * dt / dy * (i!(vn[1.., 1..]) - i!(vn[..-1, 1..])))
        );

        fill!(u[0, ..] = 1.);
        fill!(u[-1, ..] = 1.);
        fill!(u[.., 0] = 1.);
        fill!(u[.., -1] = 1.);
        fill!(v[0, ..] = 1.);
        fill!(v[-1, ..] = 1.);
        fill!(v[.., 0] = 1.);
        fill!(v[.., -1] = 1.);
    }
}
