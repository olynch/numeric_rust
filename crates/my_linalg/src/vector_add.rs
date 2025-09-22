use nalgebra::*;

pub fn for_loop_add(a: DVectorView<f64>, b: DVectorView<f64>, mut c: DVectorViewMut<f64>) {
    for i in 0..a.len() {
        c[i] = a[i] + b[i];
    }
}

pub fn for_loop_add_zip(a: DVectorView<f64>, b: DVectorView<f64>, mut c: DVectorViewMut<f64>) {
    for ((a, b), c) in a.iter().zip(b.iter()).zip(c.iter_mut()) {
        *c = a + b;
    }
}

#[inline(never)]
pub fn slice_add_to(a: &[f64], b: &[f64], c: &mut [f64]) {
    let n = c.len();
    assert_eq!(a.len(), n);
    assert_eq!(b.len(), n);
    for ((a, b), c) in a.iter().zip(b.iter()).zip(c.iter_mut()) {
        *c = a + b;
    }
}

pub fn for_loop_add_slice_zip(
    a: DVectorView<f64>,
    b: DVectorView<f64>,
    mut c: DVectorViewMut<f64>,
) {
    let a = a.as_slice();
    let b = b.as_slice();
    let c = c.as_mut_slice();
    slice_add_to(a, b, c);
}

pub fn nalgebra_add(a: DVectorView<f64>, b: DVectorView<f64>, mut c: DVectorViewMut<f64>) {
    a.add_to(&b, &mut c);
}
