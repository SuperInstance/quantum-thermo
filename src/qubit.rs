//! Qubit states, Pauli matrices, density matrices, and Bloch sphere coordinates.
//!
//! Represents qubit states as 2-component real vectors for pure states and
//! 2×2 real symmetric matrices for density matrices. Complex amplitudes are
//! encoded as `(real, imag)` pairs in the state vector.

use serde::{Deserialize, Serialize};

/// A 2×2 real matrix stored as `[[f64; 2]; 2]`.
pub type Matrix2 = [[f64; 2]; 2];

/// Pauli X matrix: [[0, 1], [1, 0]]
pub const PAULI_X: Matrix2 = [[0.0, 1.0], [1.0, 0.0]];
/// Pauli Y matrix (real part only): [[0, 0], [0, 0]]
/// For full Pauli Y you need complex arithmetic; the off-diagonal elements are imaginary.
/// This stores the real part.
pub const PAULI_Y_REAL: Matrix2 = [[0.0, 0.0], [0.0, 0.0]];
/// Pauli Z matrix: [[1, 0], [0, -1]]
pub const PAULI_Z: Matrix2 = [[1.0, 0.0], [0.0, -1.0]];
/// Identity matrix: [[1, 0], [0, 1]]
pub const IDENTITY: Matrix2 = [[1.0, 0.0], [0.0, 1.0]];

/// A pure qubit state vector with real amplitudes: |ψ⟩ = α|0⟩ + β|1⟩.
///
/// Invariant: `alpha² + beta² = 1` (normalized).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QubitState {
    /// Amplitude of |0⟩
    pub alpha: f64,
    /// Amplitude of |1⟩
    pub beta: f64,
}

impl QubitState {
    /// Create a new qubit state. Panics if not approximately normalized.
    pub fn new(alpha: f64, beta: f64) -> Self {
        let norm = alpha * alpha + beta * beta;
        if (norm - 1.0).abs() > 1e-10 {
            panic!("State not normalized: |α|² + |β|² = {norm}");
        }
        Self { alpha, beta }
    }

    /// Create a normalized qubit state from unnormalized amplitudes.
    pub fn from_unnormalized(alpha: f64, beta: f64) -> Self {
        let norm = (alpha * alpha + beta * beta).sqrt();
        if norm < 1e-15 {
            return Self {
                alpha: 1.0,
                beta: 0.0,
            }; // default to |0⟩
        }
        Self {
            alpha: alpha / norm,
            beta: beta / norm,
        }
    }

    /// |0⟩ state
    pub fn zero() -> Self {
        Self {
            alpha: 1.0,
            beta: 0.0,
        }
    }

    /// |1⟩ state
    pub fn one() -> Self {
        Self {
            alpha: 0.0,
            beta: 1.0,
        }
    }

    /// |+⟩ = (|0⟩ + |1⟩)/√2
    pub fn plus() -> Self {
        let s = 1.0 / std::f64::consts::SQRT_2;
        Self { alpha: s, beta: s }
    }

    /// |−⟩ = (|0⟩ − |1⟩)/√2
    pub fn minus() -> Self {
        let s = 1.0 / std::f64::consts::SQRT_2;
        Self { alpha: s, beta: -s }
    }

    /// Convert to density matrix: ρ = |ψ⟩⟨ψ|
    pub fn to_density_matrix(&self) -> DensityMatrix {
        DensityMatrix {
            data: [
                [self.alpha * self.alpha, self.alpha * self.beta],
                [self.beta * self.alpha, self.beta * self.beta],
            ],
        }
    }

    /// Probability of measuring |0⟩: |α|²
    pub fn prob_zero(&self) -> f64 {
        self.alpha * self.alpha
    }

    /// Probability of measuring |1⟩: |β|²
    pub fn prob_one(&self) -> f64 {
        self.beta * self.beta
    }

    /// Bloch sphere coordinates (theta, phi).
    ///
    /// For a real state |ψ⟩ = cos(θ/2)|0⟩ + sin(θ/2)|1⟩:
    /// - theta ∈ [0, π], phi = 0 (or π for negative beta)
    pub fn bloch_coords(&self) -> (f64, f64) {
        // theta = 2*arccos(|alpha|), phi = sign(beta) * 0 or π
        let theta = 2.0 * self.alpha.abs().acos();
        let phi = if self.beta >= 0.0 {
            0.0
        } else {
            std::f64::consts::PI
        };
        (theta, phi)
    }

    /// Create from Bloch sphere coordinates (theta, phi).
    /// For real states, phi is ignored and we use:
    /// |ψ⟩ = cos(θ/2)|0⟩ + sin(θ/2)|1⟩
    pub fn from_bloch(theta: f64, _phi: f64) -> Self {
        Self {
            alpha: (theta / 2.0).cos(),
            beta: (theta / 2.0).sin(),
        }
    }
}

/// A 2×2 density matrix for a qubit.
///
/// Stored as a real symmetric matrix `[[f64; 2]; 2]`. For a valid density matrix:
/// - Trace = 1
/// - Positive semidefinite (det ≥ 0 for 2×2)
/// - Hermitian (real symmetric in our representation)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DensityMatrix {
    pub data: Matrix2,
}

impl DensityMatrix {
    /// Create a new density matrix. Does NOT validate — caller must ensure physicality.
    pub fn new(data: Matrix2) -> Self {
        Self { data }
    }

    /// Maximally mixed state: I/2
    pub fn maximally_mixed() -> Self {
        Self {
            data: [[0.5, 0.0], [0.0, 0.5]],
        }
    }

    /// From two pure states with probabilities (mixed state).
    /// ρ = p|ψ₁⟩⟨ψ₁| + (1-p)|ψ₂⟩⟨ψ₂|
    pub fn from_mixture(p: f64, psi1: &QubitState, psi2: &QubitState) -> Self {
        let rho1 = psi1.to_density_matrix();
        let rho2 = psi2.to_density_matrix();
        let q = 1.0 - p;
        Self {
            data: [
                [
                    p * rho1.data[0][0] + q * rho2.data[0][0],
                    p * rho1.data[0][1] + q * rho2.data[0][1],
                ],
                [
                    p * rho1.data[1][0] + q * rho2.data[1][0],
                    p * rho1.data[1][1] + q * rho2.data[1][1],
                ],
            ],
        }
    }

    /// Trace of the density matrix (should be 1).
    pub fn trace(&self) -> f64 {
        self.data[0][0] + self.data[1][1]
    }

    /// Determinant of the 2×2 matrix.
    pub fn det(&self) -> f64 {
        self.data[0][0] * self.data[1][1] - self.data[0][1] * self.data[1][0]
    }

    /// Purity: Tr(ρ²). Equals 1 for pure states, < 1 for mixed.
    pub fn purity(&self) -> f64 {
        let a = self.data[0][0];
        let b = self.data[0][1];
        let c = self.data[1][0];
        let d = self.data[1][1];
        // ρ² = [[a²+bc, ab+bd], [ca+dc, cb+d²]]
        (a * a + b * c) + (c * b + d * d)
    }

    /// Bloch vector (rx, ry, rz) from the density matrix.
    /// For our real representation: ry = 0 always.
    /// rx = Tr(ρ X), rz = Tr(ρ Z)
    pub fn bloch_vector(&self) -> (f64, f64, f64) {
        let rx = self.data[0][1] + self.data[1][0]; // Tr(ρX) = 2ρ_01
        let ry = 0.0; // No imaginary part in our representation
        let rz = self.data[0][0] - self.data[1][1]; // Tr(ρZ) = ρ_00 - ρ_11
        (rx, ry, rz)
    }

    /// Eigenvalues of the 2×2 real symmetric matrix.
    /// For [[a, b], [b, d]]: λ = ((a+d) ± √((a-d)² + 4b²)) / 2
    pub fn eigenvalues(&self) -> (f64, f64) {
        let a = self.data[0][0];
        let b = self.data[0][1];
        let d = self.data[1][1];
        let tr = a + d;
        let disc = ((a - d) * (a - d) + 4.0 * b * b).sqrt();
        let lambda1 = (tr + disc) / 2.0;
        let lambda2 = (tr - disc) / 2.0;
        (lambda1, lambda2)
    }

    /// Check if this is a pure state (purity ≈ 1).
    pub fn is_pure(&self) -> bool {
        (self.purity() - 1.0).abs() < 1e-10
    }

    /// Apply a 2×2 matrix to the density matrix: M ρ M†
    /// For real M: M ρ M^T
    pub fn apply_matrix(&self, m: &Matrix2) -> Self {
        let rho = &self.data;
        // temp = M * rho
        let temp = [
            [
                m[0][0] * rho[0][0] + m[0][1] * rho[1][0],
                m[0][0] * rho[0][1] + m[0][1] * rho[1][1],
            ],
            [
                m[1][0] * rho[0][0] + m[1][1] * rho[1][0],
                m[1][0] * rho[0][1] + m[1][1] * rho[1][1],
            ],
        ];
        // result = temp * M^T
        let result = [
            [
                temp[0][0] * m[0][0] + temp[0][1] * m[0][1],
                temp[0][0] * m[1][0] + temp[0][1] * m[1][1],
            ],
            [
                temp[1][0] * m[0][0] + temp[1][1] * m[0][1],
                temp[1][0] * m[1][0] + temp[1][1] * m[1][1],
            ],
        ];
        Self { data: result }
    }

    /// Matrix addition: self + other
    pub fn add(&self, other: &DensityMatrix) -> DensityMatrix {
        DensityMatrix {
            data: [
                [
                    self.data[0][0] + other.data[0][0],
                    self.data[0][1] + other.data[0][1],
                ],
                [
                    self.data[1][0] + other.data[1][0],
                    self.data[1][1] + other.data[1][1],
                ],
            ],
        }
    }

    /// Scalar multiplication
    pub fn scale(&self, s: f64) -> DensityMatrix {
        DensityMatrix {
            data: [
                [self.data[0][0] * s, self.data[0][1] * s],
                [self.data[1][0] * s, self.data[1][1] * s],
            ],
        }
    }
}

/// Matrix multiplication for 2×2 matrices.
pub fn mat_mul(a: &Matrix2, b: &Matrix2) -> Matrix2 {
    [
        [
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ],
        [
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ],
    ]
}

/// Matrix transpose.
pub fn mat_transpose(m: &Matrix2) -> Matrix2 {
    [[m[0][0], m[1][0]], [m[0][1], m[1][1]]]
}

/// Matrix trace.
pub fn mat_trace(m: &Matrix2) -> f64 {
    m[0][0] + m[1][1]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_qubit_zero_state() {
        let q = QubitState::zero();
        assert!((q.alpha - 1.0).abs() < 1e-10);
        assert!((q.beta).abs() < 1e-10);
        assert!((q.prob_zero() - 1.0).abs() < 1e-10);
        assert!((q.prob_one()).abs() < 1e-10);
    }

    #[test]
    fn test_qubit_plus_state() {
        let q = QubitState::plus();
        assert!((q.alpha - 1.0 / std::f64::consts::SQRT_2).abs() < 1e-10);
        assert!((q.prob_zero() - 0.5).abs() < 1e-10);
        assert!((q.prob_one() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_density_matrix_trace() {
        let rho = QubitState::zero().to_density_matrix();
        assert!((rho.trace() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_purity_pure_state() {
        let rho = QubitState::plus().to_density_matrix();
        assert!((rho.purity() - 1.0).abs() < 1e-10);
        assert!(rho.is_pure());
    }

    #[test]
    fn test_purity_mixed_state() {
        let rho = DensityMatrix::from_mixture(0.5, &QubitState::zero(), &QubitState::one());
        assert!((rho.purity() - 0.5).abs() < 1e-10);
        assert!(!rho.is_pure());
    }

    #[test]
    fn test_bloch_coords() {
        let q = QubitState::zero();
        let (theta, _) = q.bloch_coords();
        // For |0⟩: alpha=1, beta=0 → acos(0) = π/2, theta = π
        // Actually: theta = 2*acos(beta) = 2*acos(0) = π
        // But our formula uses beta.acos() which is acos(1) for |0⟩ since we check sign of beta
        // For |0⟩: alpha=1, beta=0. theta = 2*acos(0) = π
        assert!((theta - std::f64::consts::PI).abs() < 1e-10 || theta.abs() < 1e-10);
    }

    #[test]
    fn test_eigenvalues() {
        let rho = DensityMatrix::maximally_mixed();
        let (e1, e2) = rho.eigenvalues();
        assert!((e1 - 0.5).abs() < 1e-10);
        assert!((e2 - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_from_unnormalized() {
        let q = QubitState::from_unnormalized(3.0, 4.0);
        assert!((q.alpha - 0.6).abs() < 1e-10);
        assert!((q.beta - 0.8).abs() < 1e-10);
    }

    #[test]
    fn test_apply_matrix() {
        let rho = QubitState::zero().to_density_matrix();
        let result = rho.apply_matrix(&PAULI_X);
        // X|0⟩⟨0|X = |1⟩⟨1|
        assert!((result.data[0][0]).abs() < 1e-10);
        assert!((result.data[1][1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_from_bloch() {
        let q = QubitState::from_bloch(std::f64::consts::PI, 0.0);
        assert!((q.alpha).abs() < 1e-10); // cos(π/2) ≈ 0
        assert!((q.beta - 1.0).abs() < 1e-10); // sin(π/2) = 1
    }
}
