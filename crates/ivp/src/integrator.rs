use nalgebra::*;
use num_traits::Float;

use crate::{
    adaptive_strategy::AdaptiveStrategy, step_algorithm::StepAlgorithm, system::OdeSystem,
};

#[derive(Clone, Copy)]
pub struct TSpan<F> {
    pub start: F,
    pub end: F,
}

impl<F> TSpan<F> {
    pub fn new(start: F, end: F) -> Self {
        Self { start, end }
    }
}

pub struct Integrator<
    'a,
    F: Scalar + Float,
    Sys: OdeSystem<F>,
    Step: StepAlgorithm<F>,
    AS: AdaptiveStrategy<F, Step::ErrorEstimate>,
> {
    pub(crate) sys: &'a Sys,
    pub(crate) step_algorithm: &'a Step,
    pub(crate) adaptive_strategy: &'a AS,
    pub(crate) dt: F,
    pub(crate) cache: Step::Cache,
    pub(crate) k: usize,
    pub(crate) ts: Vec<F>,
    pub(crate) ys: Vec<DVector<F>>,
    pub(crate) interpolants: Vec<Step::Interpolant>,
}

impl<
    'a,
    F: Float + Scalar,
    Sys: OdeSystem<F>,
    Step: StepAlgorithm<F>,
    AS: AdaptiveStrategy<F, Step::ErrorEstimate>,
> Integrator<'a, F, Sys, Step, AS>
{
    pub fn new(
        sys: &'a Sys,
        step_algorithm: &'a Step,
        adaptive_strategy: &'a AS,
        tspan: TSpan<F>,
        y0: DVector<F>,
    ) -> Self {
        let cache = step_algorithm.init_cache(sys);
        let dt = adaptive_strategy.init_dt();
        let t0 = tspan.start;
        Self {
            sys,
            step_algorithm,
            adaptive_strategy,
            dt,
            cache,
            k: 0,
            ts: vec![t0],
            ys: vec![y0],
            interpolants: vec![],
        }
    }

    pub fn step(&mut self) {
        let y0 = self.ys[self.k].as_view();
        let mut y1 = DVector::zeros(self.sys.dimension());
        let t = self.ts[self.k];
        let (interpolant, error) =
            self.step_algorithm
                .step(&mut self.cache, self.sys, y1.as_view_mut(), y0, t, self.dt);
        match self.adaptive_strategy.try_accept(self.dt, error, y0) {
            Ok(dt) => {
                self.ts.push(t + self.dt);
                self.dt = dt;
                self.k += 1;
                self.ys.push(y1);
                self.interpolants.push(interpolant);
            }
            Err(dt) => {
                self.dt = dt;
                self.step();
            }
        }
    }
}
