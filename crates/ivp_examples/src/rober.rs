use ivp::*;
use nalgebra::*;

#[derive(Clone, Copy)]
pub struct Rober<F> {
    k1: F,
    k2: F,
    k3: F,
}

impl OdeSystem<f64> for Rober<f64> {
    fn dimension(&self) -> usize {
        3
    }

    fn labels(&self) -> Vec<String> {
        vec!["y1".to_string(), "y2".to_string(), "y3".to_string()]
    }

    fn vfield(&self, mut du: DVectorViewMut<f64>, u: DVectorView<f64>, _t: f64) {
        let Rober { k1, k2, k3 } = *self;
        let (y1, y2, y3) = (u[0], u[1], u[2]);
        du[0] = -k1 * y1 + k3 * y2 * y3;
        du[1] = k1 * y1 - k3 * y2 * y3 - k2 * y2 * y2;
        du[2] = y1 + y2 + y3 - 1.;
    }

    fn mass_matrix(&self) -> DMatrix<f64> {
        dmatrix![
            1., 0., 0.;
            0., 1., 0.;
            0., 0., 0.
        ]
    }

    fn jacobian(&self, mut out: DMatrixViewMut<f64>, u: DVectorView<f64>, _t: f64) {
        let Rober { k1, k2, k3 } = *self;
        let (_y1, y2, y3) = (u[0], u[1], u[2]);
        out[(0, 0)] = -k1;
        out[(0, 1)] = k3 * y3;
        out[(0, 2)] = k3 * y2;
        out[(1, 0)] = k1;
        out[(1, 1)] = -k3 * y3 - 2.0 * k2 * y2;
        out[(1, 2)] = -k3 * y2;
        out[(2, 0)] = 1.0;
        out[(2, 1)] = 1.0;
        out[(2, 2)] = 1.0;
    }
}

const EXAMPLE_SYS: Rober<f64> = Rober {
    k1: 0.04,
    k2: 3e7,
    k3: 1e4,
};

pub fn create_prob() -> OdeProblem<f64, Rober<f64>> {
    OdeProblem::new(EXAMPLE_SYS, dvector![1.0, 0.0, 0.0], TSpan::new(1e-5, 1e5))
}
