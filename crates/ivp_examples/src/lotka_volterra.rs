use std::fmt::Debug;

use ivp::*;
use nalgebra::*;
use num_traits::Float;

#[derive(Clone, Copy)]
pub struct LotkaVolterra<F> {
    alpha: F,
    beta: F,
    gamma: F,
    delta: F,
}

impl<F: Float + Debug + 'static> OdeSystem<F> for LotkaVolterra<F> {
    fn dimension(&self) -> usize {
        2
    }

    fn labels(&self) -> Vec<String> {
        vec!["sheep".to_string(), "wolves".to_string()]
    }

    fn vfield(&self, mut out: DVectorViewMut<F>, y: DVectorView<F>, _t: F) {
        let LotkaVolterra {
            alpha,
            beta,
            gamma,
            delta,
        } = *self;
        let (x, y) = (y[0], y[1]);
        out[0] = alpha * x - beta * x * y;
        out[1] = delta * x * y - gamma * y;
    }
}

const EXAMPLE_SYS: LotkaVolterra<f64> = LotkaVolterra {
    alpha: 1.5,
    beta: 1.0,
    gamma: 3.0,
    delta: 1.0,
};

pub fn create_prob() -> OdeProblem<f64, LotkaVolterra<f64>> {
    OdeProblem::new(EXAMPLE_SYS, dvector![1.0, 1.0], TSpan::new(0.0, 10.0))
}
