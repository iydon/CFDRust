use cfd_rust::lesson;

#[test]
fn x1() {
    let (mut u, nt, dx, dt, c) = lesson::x1::default();
    lesson::x1::solve(&mut u, nt, dx, dt, c);
}

#[test]
fn x2() {
    let (mut u, nt, dx, dt) = lesson::x2::default();
    lesson::x2::solve(&mut u, nt, dx, dt);
}

#[test]
fn x3() {
    let (mut u, nt, dx, dt, c) = lesson::x3::default();
    lesson::x3::solve(&mut u, nt, dx, dt, c);
}

#[test]
fn x4() {
    let (mut u, nt, dx, dt, nu) = lesson::x4::default();
    lesson::x4::solve(&mut u, nt, dx, dt, nu);
}

#[test]
fn x5() {
    let (mut u, nt, dx, dt, nu) = lesson::x5::default();
    lesson::x5::solve(&mut u, nt, dx, dt, nu);
}

#[test]
fn x7() {
    let (mut u, nt, dx, dy, dt, c) = lesson::x7::default();
    lesson::x7::solve(&mut u, nt, dx, dy, dt, c);
}

#[test]
fn x8() {
    let (mut u, mut v, nt, dx, dy, dt) = lesson::x8::default();
    lesson::x8::solve(&mut u, &mut v, nt, dx, dy, dt);
}

#[test]
fn x9() {
    let (mut u, nt, dx, dy, dt, nu) = lesson::x9::default();
    lesson::x9::solve(&mut u, nt, dx, dy, dt, nu);
}

#[test]
fn xa() {
    let (mut u, mut v, nt, dx, dy, dt, nu) = lesson::xa::default();
    lesson::xa::solve(&mut u, &mut v, nt, dx, dy, dt, nu);
}

#[test]
fn xc() {
    let (mut p, y, dx, dy, eps) = lesson::xc::default();
    lesson::xc::solve(&mut p, y, dx, dy, eps);
}

#[test]
fn xd() {
    let (mut p, b, nt, dx, dy) = lesson::xd::default();
    lesson::xd::solve(&mut p, b, nt, dx, dy);
}

#[test]
fn xe() {
    let (mut u, mut v, mut p, nt, nit, dx, dy, dt, rho, nu) = lesson::xe::default();
    lesson::xe::solve(&mut u, &mut v, &mut p, nt, nit, dx, dy, dt, rho, nu);
}

#[test]
fn xf() {
    let (mut u, mut v, mut p, nit, dx, dy, dt, rho, nu, f, eps) = lesson::xf::default();
    lesson::xf::solve(&mut u, &mut v, &mut p, nit, dx, dy, dt, rho, nu, f, eps);
}
