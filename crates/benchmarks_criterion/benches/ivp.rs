use criterion::*;
use ivp::*;
use ivp_examples::*;
use nalgebra::*;
use num_traits::Float;

fn bench_ivp<
    F: Float + Scalar + 'static,
    S: OdeSystem<F>,
    SA: StepAlgorithm<F>,
    AS: AdaptiveStrategy<F, SA::ErrorEstimate>,
>(
    name: &str,
    c: &mut Criterion,
    prob: OdeProblem<F, S>,
    step_algorithm: SA,
    adaptive_strategy: AS,
) {
    c.bench_function(name, |b| {
        b.iter(|| {
            prob.solve(&step_algorithm, &adaptive_strategy);
        })
    });
}

fn ivp(c: &mut Criterion) {
    bench_ivp(
        "lotka-volterra tsit5",
        c,
        lotka_volterra::create_prob(),
        Tsit5,
        AdaptiveStep::new(0.1, 1e-3, 1e-6, 4),
    );

    // bench_ivp(
    //     "lotka-volterra rosenbrock23",
    //     c,
    //     lotka_volterra::create_prob(),
    //     Rosenbrock23,
    //     AdaptiveStep::new(0.1, 1e-3, 1e-6, 4),
    // );

    bench_ivp(
        "pleiades tsit5",
        c,
        pleiades::create_prob(),
        Tsit5,
        AdaptiveStep::new(0.1, 1e-3, 1e-6, 4),
    );
}

criterion_group!(benches, ivp);
criterion_main!(benches);
