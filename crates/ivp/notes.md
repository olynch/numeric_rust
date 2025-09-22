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

# Optimization log

## Tsit5

### Implemented optimizations

- Store the `ks` as `Matrix<F, U7, Dyn, ...>` instead of `[DVector<F>; 7]`. This gave 60% speedup on lotka-volterra, 20% speedup on pleiades.
- Use [crunchy](https://github.com/eira-fransham/crunchy) to manually unroll the stage loop. This didn't make a difference for lotka-volterra, but gave a 15% speedup on pleiades. This implies to me that LLVM already decided to unroll the loop for lotka-volterra, but didn't unroll the loop for pleiades.
- Change the order of iteration for the computation of y1hat to put `s` in the outside loop. This seems to at best give a marginal improvement.
- Use [mimalloc](https://github.com/microsoft/mimalloc). This gives a 10-30% speedup, essentially for free. That is pretty fantastic; I'm going to use mimalloc in all of my Rust programs now. This does imply that using a bump allocator would probably be even faster. Oddly enough, mimalloc v2 seems faster than v3.
