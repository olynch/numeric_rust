use ivp::*;
use nalgebra::*;
use num_traits::Float;

pub struct Pleaides {
    nplanets: usize,
}

impl<F: Float + Scalar> OdeSystem<F> for Pleaides {
    fn dimension(&self) -> usize {
        self.nplanets * 4
    }

    fn labels(&self) -> Vec<String> {
        let n = self.nplanets;
        (0..n)
            .map(|i| format!("x{i}"))
            .chain((0..n).map(|i| format!("y{i}")))
            .chain((0..n).map(|i| format!("vx{i}")))
            .chain((0..n).map(|i| format!("vy{i}")))
            .collect()
    }

    fn vfield(&self, mut out: DVectorViewMut<F>, u: DVectorView<F>, _t: F) {
        let n = self.nplanets;
        let x = &u.as_slice()[0..n];
        let y = &u.as_slice()[n..2 * n];

        out.view_mut((0, 0), (2 * n, 1))
            .copy_from(&u.view((2 * n, 0), (2 * n, 1)));

        for i in 0..n {
            for j in 0..n {
                if i != j {
                    let dx = x[j] - x[i];
                    let dy = y[j] - y[i];
                    let r = (dx * dx + dy * dy).sqrt();
                    let r3 = r * r * r;
                    out[2 * n + i] = out[2 * n + i] + F::from(j).unwrap() * dx / r3;
                    out[3 * n + i] = out[3 * n + i] + F::from(j).unwrap() * dy / r3;
                }
            }
        }
    }
}

const EXAMPLE_SYS: Pleaides = Pleaides { nplanets: 7 };

pub fn create_prob() -> OdeProblem<f64, Pleaides> {
    OdeProblem::new(
        EXAMPLE_SYS,
        dvector![
            3.0, 3.0, -1.0, -3.0, 2.0, -2.0, 2.0, 3.0, -3.0, 2.0, 0., 0., -4.0, 4.0, 0., 0., 0.,
            0., 0., 1.75, -1.5, 0., 0., 0., -1.25, 1., 0., 0.
        ],
        TSpan::new(0.0, 3.0),
    )
}
