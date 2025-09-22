use gungraun::*;
use my_linalg::vector_add::*;
use nalgebra::DVector;

fn generate_vectors() -> (DVector<f64>, DVector<f64>, DVector<f64>) {
    let size = 512;
    let a = DVector::<f64>::new_random(size);
    let b = DVector::<f64>::new_random(size);
    let c = DVector::<f64>::zeros(size);
    (a, b, c)
}

#[library_benchmark]
#[bench::first(args = (), setup = generate_vectors)]
fn for_loop_add_bench(
    mut abc: (DVector<f64>, DVector<f64>, DVector<f64>),
) -> (DVector<f64>, DVector<f64>, DVector<f64>) {
    for_loop_add(abc.0.as_view(), abc.1.as_view(), abc.2.as_view_mut());
    abc
}

#[library_benchmark]
#[bench::first(args = (), setup = generate_vectors)]
fn for_loop_add_slice_zip_bench(
    mut abc: (DVector<f64>, DVector<f64>, DVector<f64>),
) -> (DVector<f64>, DVector<f64>, DVector<f64>) {
    for_loop_add_slice_zip(abc.0.as_view(), abc.1.as_view(), abc.2.as_view_mut());
    abc
}

library_benchmark_group!(name = vector_add; benchmarks = for_loop_add_bench, for_loop_add_slice_zip_bench);
main!(library_benchmark_groups = vector_add);
