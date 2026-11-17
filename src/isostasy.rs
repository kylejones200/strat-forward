//! Deflection of a thin elastic plate under a load, on the model's 1D grid.
//!
//! Solves `(D B + diag(k)) w = q`, where `B` is the bending operator, `k` the restoring force per unit
//! deflection (Δρ g, per cell) and `q` the load per unit area (Pa). `B` comes from the plate's bending
//! energy `D/2 ∫ w''² dx` with curvature taken at interior nodes only, which leaves both ends free (no
//! moment, no shear). It is symmetric and annihilates any straight line, so a uniform load gives a uniform
//! deflection and the solution balances total force and moment exactly:
//! `Σ k w = Σ q` and `Σ x k w = Σ x q`. Interior rows are the standard `[1, −4, 6, −4, 1] / dx⁴`.

/// Standard gravity (m/s²).
pub const G: f64 = 9.81;

/// Flexural parameter α = (4 D / (Δρ g))^¼ (m): the length over which a line load's deflection decays.
pub fn flexural_parameter(rigidity: f64, restoring: f64) -> f64 {
    (4.0 * rigidity / restoring).powf(0.25)
}

/// Deflection `w` (m, positive down) of a plate of rigidity `rigidity` (N·m) on grid spacing `dx` (m) under
/// load `load` (Pa) with restoring force `restoring` (Pa/m) per cell. `restoring` must be positive.
pub fn deflection(dx: f64, rigidity: f64, restoring: &[f64], load: &[f64]) -> Vec<f64> {
    let n = load.len();
    assert_eq!(restoring.len(), n);
    // Symmetric pentadiagonal matrix: main diagonal, first and second super-diagonals.
    let mut a0 = restoring.to_vec();
    let mut a1 = vec![0.0; n];
    let mut a2 = vec![0.0; n];
    if n >= 3 && rigidity > 0.0 {
        let s = rigidity / dx.powi(4);
        // Each interior node r contributes the outer product of the curvature stencil (1, −2, 1) at r−1..r+1.
        const C: [f64; 3] = [1.0, -2.0, 1.0];
        for r in 1..n - 1 {
            for p in 0..3 {
                let i = r - 1 + p;
                a0[i] += s * C[p] * C[p];
                if p + 1 < 3 {
                    a1[i] += s * C[p] * C[p + 1];
                }
                if p + 2 < 3 {
                    a2[i] += s * C[p] * C[p + 2];
                }
            }
        }
    }
    // LDLᵀ without pivoting: the matrix is symmetric positive definite (bending is semi-definite, k > 0).
    let (mut d, mut l1, mut l2) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    for i in 0..n {
        if i >= 2 {
            l2[i] = a2[i - 2] / d[i - 2];
        }
        if i >= 1 {
            let back = if i >= 2 {
                l2[i] * l1[i - 1] * d[i - 2]
            } else {
                0.0
            };
            l1[i] = (a1[i - 1] - back) / d[i - 1];
        }
        let mut di = a0[i];
        if i >= 1 {
            di -= l1[i] * l1[i] * d[i - 1];
        }
        if i >= 2 {
            di -= l2[i] * l2[i] * d[i - 2];
        }
        d[i] = di;
    }
    let mut w = load.to_vec();
    for i in 0..n {
        if i >= 1 {
            w[i] -= l1[i] * w[i - 1];
        }
        if i >= 2 {
            w[i] -= l2[i] * w[i - 2];
        }
    }
    for (wi, di) in w.iter_mut().zip(&d) {
        *wi /= di;
    }
    for i in (0..n).rev() {
        if i + 1 < n {
            w[i] -= l1[i + 1] * w[i + 1];
        }
        if i + 2 < n {
            w[i] -= l2[i + 2] * w[i + 2];
        }
    }
    w
}
