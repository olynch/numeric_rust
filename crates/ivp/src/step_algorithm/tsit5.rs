use crate::OdeSystem;

use super::StepAlgorithm;
use crunchy::unroll;
use nalgebra::*;
use num_traits::Float;

const C: [f64; 7] = [0.0, 0.161, 0.327, 0.9, 0.9800255409045097, 1., 1.];

const B: [f64; 7] = [
    0.09646076681806523,
    0.01,
    0.4798896504144996,
    1.379008574103742,
    -3.290069515436081,
    2.324710524099774,
    0.0,
];

const B_TILDE: [f64; 7] = [
    0.001780011052226,
    0.000816434459657,
    -0.007880878010262,
    0.144711007173263,
    -0.582357165452555,
    0.458082105929187,
    1. / 66.,
];

#[rustfmt::skip]
const A: [[f64; 7]; 7] = {
    let a: [[f64; 7]; 6] = [
        [0., 0., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0.],
        [0., 0.3354806554923570, 0., 0., 0., 0., 0.],
        [0., -6.359448489975075, 4.362295432869581, 0., 0., 0., 0.],
        [0., -11.74888356406283, 7.495539342889836, -0.09249506636175525, 0., 0., 0.],
        [0., -12.92096931784711, 8.159367898576159, -0.07158497328140100, -0.02826905039406838, 0., 0.],
    ];
    [
        [0.,                                           0.,      0.,      0.,      0.,      0., 0.],
        [C[1],                                         0.,      0.,      0.,      0.,      0., 0.],
        [C[2] - a[2][1],                               a[2][1], 0.,      0.,      0.,      0., 0.],
        [C[3] - a[3][1] - a[3][2],                     a[3][1], a[3][2], 0.,      0.,      0., 0.],
        [C[4] - a[4][1] - a[4][2] - a[4][3],           a[4][1], a[4][2], a[4][3], 0.,      0., 0.],
        [C[5] - a[5][1] - a[5][2] - a[5][3] - a[5][4], a[5][1], a[5][2], a[5][3], a[5][4], 0., 0.],
        B
    ]
};

pub struct Tsit5Cache<F> {
    c: [F; 7],
    btilde: [F; 7],
    a: [[F; 7]; 7],
}

pub struct Tsit5;

impl<F: Float + Scalar + ComplexField<RealField = F>> StepAlgorithm<F> for Tsit5 {
    type Cache = Tsit5Cache<F>;
    type Interpolant = Matrix<F, Dyn, U7, VecStorage<F, Dyn, U7>>;
    type ErrorEstimate = F;

    fn init_cache<S: OdeSystem<F>>(&self, _sys: &S) -> Self::Cache {
        Tsit5Cache {
            c: C.map(|x| F::from(x).unwrap()),
            btilde: B_TILDE.map(|x| F::from(x).unwrap()),
            a: A.map(|r| r.map(|x| F::from(x).unwrap())),
        }
    }

    fn interpolate(
        &self,
        y0: DVectorView<F>,
        _y1: DVectorView<F>,
        ks: &Self::Interpolant,
        dt: F,
        t: F,
    ) -> DVector<F> {
        let f = |x| F::from(x).unwrap();
        let t2 = t * t;
        let btilde_t = [
            f(-1.0530884977290216)
                * t
                * (t - f(1.3299890189751412))
                * (t2 - f(1.4364028541716351) * t + f(0.7139816917074209)),
            f(0.1017) * t2 * (t2 - f(2.1966568338249754) * t + f(1.2949852507374631)),
            f(2.490627285651252793)
                * t2
                * (t2 - f(2.38535645472061657) * t + f(1.57803468208092486)),
            f(-16.54810288924490272)
                * (t - f(1.21712927295533244))
                * (t - f(0.61620406037800089))
                * t2,
            f(47.37952196281928122)
                * (t - f(1.203071208372362603))
                * (t - f(0.658047292653547382))
                * t2,
            f(-34.87065786149660974) * (t - f(1.2)) * (t - f(0.666666666666666667)) * t2,
            f(2.5) * (t - f(1.)) * (t - f(0.6)) * t2,
        ];
        let mut y = y0.into_owned();
        for i in 0..y0.len() {
            y[i] = y[i]
                + dt * (0..7)
                    .map(|s| btilde_t[s] * ks[(i, s)])
                    .fold(F::zero(), |x, y| x + y);
        }
        y
    }

    fn step<S: crate::OdeSystem<F>>(
        &self,
        cache: &mut Self::Cache,
        system: &S,
        mut y1: DVectorViewMut<F>,
        y0: DVectorView<F>,
        t: F,
        dt: F,
    ) -> (Self::Interpolant, Self::ErrorEstimate) {
        let n = system.dimension();
        let mut ks = Self::Interpolant::zeros(n);
        y1.copy_from(&y0);
        unroll! { for s in 0..7 {
            system.vfield(ks.column_mut(s), y1.as_view(), t + cache.c[s] * dt);
            for i in 0..n {
                let mut dy = F::zero();
                for j in 0..s {
                    dy += cache.a[s][j] * ks[(i, j)];
                }
                y1[i] = y0[i] + dt * dy;
            }
        } }
        let mut y1hat = y0.into_owned();
        unroll! { for s in 0..7 {
            let btilde_dt = cache.btilde[s] * dt;
            for i in 0..n {
                y1hat[i] = y1hat[i] + btilde_dt * ks[(i, s)];
            }
        }}
        let error = EuclideanNorm.metric_distance(&y1, &y1hat);
        (ks, error)
    }
}
