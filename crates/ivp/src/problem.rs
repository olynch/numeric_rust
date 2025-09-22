use std::fmt::Write as _;
use std::fmt::{Debug, Display};
use std::io::Write as _;
use std::process::Command;
use std::process::Stdio;

use crate::{AdaptiveStrategy, Integrator, OdeSystem, StepAlgorithm, TSpan};
use nalgebra::*;
use num_traits::{Float, cast};
use tempfile::NamedTempFile;

pub struct OdeProblem<F, S> {
    sys: S,
    y0: DVector<F>,
    tspan: TSpan<F>,
}

impl<F: Float + Scalar, S: OdeSystem<F>> OdeProblem<F, S> {
    pub fn new(sys: S, y0: DVector<F>, tspan: TSpan<F>) -> Self {
        Self { sys, y0, tspan }
    }

    pub fn solve<'a, SA: StepAlgorithm<F>, AS: AdaptiveStrategy<F, SA::ErrorEstimate>>(
        &'a self,
        step_algorithm: &'a SA,
        adaptive_strategy: &'a AS,
    ) -> OdeSolution<'a, F, SA> {
        let mut integrator = Integrator::new(
            &self.sys,
            step_algorithm,
            adaptive_strategy,
            self.tspan,
            self.y0.clone(),
        );
        while *integrator.ts.last().unwrap() < self.tspan.end {
            integrator.step();
        }
        OdeSolution {
            labels: self.sys.labels(),
            tspan: self.tspan,
            ts: integrator.ts,
            ys: integrator.ys,
            step_algorithm: integrator.step_algorithm,
            interpolants: integrator.interpolants,
        }
    }
}

pub struct OdeSolution<'a, F: Float + Scalar + 'static, SA: StepAlgorithm<F>> {
    labels: Vec<String>,
    tspan: TSpan<F>,
    ts: Vec<F>,
    ys: Vec<DVector<F>>,
    step_algorithm: &'a SA,
    interpolants: Vec<SA::Interpolant>,
}

impl<'a, F: Float + Scalar, SA: StepAlgorithm<F>> OdeSolution<'a, F, SA> {
    pub fn solution_at(&self, t: F) -> DVector<F> {
        if t <= self.tspan.start {
            return self.ys[0].clone();
        } else if t >= self.tspan.end {
            return self.ys.last().unwrap().clone();
        }
        let k = self.ts.partition_point(|ti| *ti < t).min(self.ts.len() - 1);
        let y0 = self.ys[k - 1].as_view();
        let y1 = self.ys[k].as_view();
        let t0 = self.ts[k - 1];
        let t1 = self.ts[k];
        let dt = t1 - t0;
        self.step_algorithm
            .interpolate(y0, y1, &self.interpolants[k - 1], dt, (t - t0) / dt)
    }

    pub fn gnuplot(&self, dt: f32)
    where
        F: Display + Debug,
    {
        let n: usize =
            cast(((self.tspan.end - self.tspan.start) / F::from(dt).unwrap()).ceil()).unwrap();
        let mut tmp = NamedTempFile::new().unwrap();
        let mut script = String::new();
        write!(&mut script, "plot ").unwrap();
        for (i, label) in self.labels.iter().enumerate() {
            if i > 0 {
                write!(&mut script, ", ").unwrap();
            }
            write!(
                &mut script,
                "'{}' using 1:{} with lines title '{}'",
                tmp.as_ref().display(),
                i + 2,
                label
            )
            .unwrap();
        }
        for k in 0..n {
            let t = F::from(k as f32 * dt).unwrap();
            write!(&mut tmp, "{t}").unwrap();
            let y = self.solution_at(t);
            for i in 0..self.labels.len() {
                write!(&mut tmp, ", {}", y[i]).unwrap();
            }
            writeln!(&mut tmp).unwrap();
        }
        let mut process = Command::new("gnuplot")
            .args(["-p", "-e", script.as_str()])
            .stdin(Stdio::piped())
            .spawn()
            .expect("failed to start plotting process");

        process.wait().expect("failed to plot");
    }

    pub fn gnuplot_log(&self, dt: F)
    where
        F: Display + Debug,
    {
        let mut tmp = NamedTempFile::new().unwrap();
        let mut script_tmp = NamedTempFile::new().unwrap();
        writeln!(&mut script_tmp, "set multiplot layout 3,1").unwrap();
        for (i, label) in self.labels.iter().enumerate() {
            writeln!(
                &mut script_tmp,
                "plot '{}' using 1:{} with lines title '{}'",
                tmp.as_ref().display(),
                i + 2,
                label
            )
            .unwrap();
        }
        writeln!(&mut script_tmp, "unset multiplot").unwrap();
        let mut t = self.tspan.start;
        let factor = F::from(10.).unwrap().powf(dt);
        while t < self.tspan.end {
            write!(&mut tmp, "{}", t.log10()).unwrap();
            // print!("{}", t.log10());
            let y = self.solution_at(t);
            for i in 0..self.labels.len() {
                write!(&mut tmp, ", {}", y[i]).unwrap();
                // print!(", {}", y[i]);
            }
            writeln!(&mut tmp).unwrap();
            // println!();
            t = t * factor;
        }
        let mut process = Command::new("gnuplot")
            .args(["-p", script_tmp.path().as_os_str().to_str().unwrap()])
            .stdin(Stdio::piped())
            .spawn()
            .expect("failed to start plotting process");

        process.wait().expect("failed to plot");
    }
}
