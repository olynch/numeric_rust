# Adaptive ODE solving and error estimates

For adaptive ODE solving, it is necessary at each step to also compute an error estimate.

Typically, this error estimate is given by comparing the current method to a lower-order method. Of course, this only makes sense when such a lower-order method exists, which is not true for Euler's method. However, for Runge-Kutta and others, this is true.

The general procedure for adaptive solving with such an error estimate is the following.

Suppose that `E` is the estimate. Then we compute a step size according to the formula

`hₙ₊₁ := hₙ(1/E)^(1/(q+1))`

where `q` is the order of the lower-order method.

# Next steps

Implement a MLIR backend for a closure-free 1ML, using type-level shapes for shape checking?

Implement a MLIR backend for an embedded DSL in Rust using dynamic shape checking. No control flow. We could use egg for this, but let's keep it simple.

Could we do this as a DSL in Lean?
