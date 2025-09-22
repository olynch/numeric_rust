use criterion::*;
use my_linalg::vector_add::*;
use nalgebra::*;

fn bench_add(
    name: &str,
    c: &mut Criterion,
    f: impl Fn(DVectorView<f64>, DVectorView<f64>, DVectorViewMut<f64>) -> (),
) {
    static KB: usize = 1024;
    let mut group = c.benchmark_group(name);
    for size in [KB, 2 * KB, 4 * KB].iter() {
        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, &size| {
            b.iter_batched_ref(
                || {
                    let a = DVector::<f64>::new_random(size);
                    let b = DVector::<f64>::new_random(size);
                    let c = DVector::<f64>::zeros(size);
                    (a, b, c)
                },
                |(a, b, c)| {
                    f(a.as_view(), b.as_view(), c.as_view_mut());
                },
                BatchSize::LargeInput,
            );
        });
    }
}

// Conclusion: looking at the source of vectorized_add, it looks like the only
// trick is to remove bounds checking and to operate on a raw slice. So if the
// loop ends up simdifying, that's purely a rust optimization thing.
//
// It seems like I can get equivalent performance to that via the
// `for_loop_add_slice_zip`, which presumably avoids bounds checking in a
// similar way.
//
// However, for_loop_add_zip, which doesn't convert the vectors into slices first
// is absolutely terrible. I think this is because the iteration logic for DVectors
// is too complicated for Rust to autovectorize? Or something similar to that.
// It gets the same 100% speedup when we switch to use native CPU features
// (which enables 4xfloat SIMD instead of 2xfloat SIMD), so it's probable using
// SIMD, but some kind of branching is slower, idk.
//
// I think the conclusion here is that the code which best expresses what is
// going on to the rust compiler (the zipped for loop over slices) is most
// efficient.
//
// It's also easiest to use this. I want to rewrite ivp to use a bump allocator
// and just work on slices.
//
// I need to figure out how to benchmark the ivp code to figure out how much
// time is spent in memory management (allocation, cache misses) vs. arithmetic.
// Is there criterion stuff for that?
//
// The next thing to benchmark is nalgebra lu solve vs. LAPACK lu solve.
pub fn addition(c: &mut Criterion) {
    bench_add("for loop add", c, for_loop_add);
    bench_add("for loop add zip", c, for_loop_add_zip);
    bench_add("for loop add slice zip", c, for_loop_add_slice_zip);
    bench_add("nalgebra add", c, nalgebra_add);
}

criterion_group!(benches, addition);
criterion_main!(benches);
