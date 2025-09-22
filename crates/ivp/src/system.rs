use nalgebra::*;
use num_traits::Float;

pub trait OdeSystem<F: Scalar + Float> {
    fn dimension(&self) -> usize;

    fn labels(&self) -> Vec<String>;

    fn vfield(&self, out: DVectorViewMut<F>, y: DVectorView<F>, t: F);

    fn mass_matrix(&self) -> DMatrix<F> {
        let n = self.dimension();
        DMatrix::identity(n, n)
    }

    fn jacobian(&self, _out: DMatrixViewMut<F>, _y: DVectorView<F>, _t: F) {
        panic!("Jacobian not implemented");
    }
}
