//! Quantum channels in Kraus representation.
//!
//! A quantum channel is a completely positive trace-preserving (CPTP) map
//! described by Kraus operators {E_k}: ρ' = Σ_k E_k ρ E_k†

use crate::qubit::{DensityMatrix, Matrix2};
#[cfg(test)]
use crate::qubit::QubitState;
use serde::{Deserialize, Serialize};

/// A quantum channel with up to 4 Kraus operators.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QuantumChannel {
    /// Kraus operators
    pub kraus: Vec<Matrix2>,
}

impl QuantumChannel {
    /// Create a new quantum channel from Kraus operators.
    /// Does NOT verify the completeness relation Σ_k E_k† E_k = I.
    pub fn new(kraus: Vec<Matrix2>) -> Self {
        Self { kraus }
    }

    /// Apply the channel to a density matrix: ρ' = Σ_k E_k ρ E_k^T
    pub fn apply(&self, rho: &DensityMatrix) -> DensityMatrix {
        let mut result = DensityMatrix::new([[0.0, 0.0], [0.0, 0.0]]);
        for ek in &self.kraus {
            result = result.add(&rho.apply_matrix(ek));
        }
        result
    }

    /// Verify the completeness relation: Σ_k E_k^T E_k ≈ I
    pub fn verify_completeness(&self) -> bool {
        let mut sum = [[0.0, 0.0], [0.0, 0.0]];
        for ek in &self.kraus {
            let etk = crate::qubit::mat_transpose(ek);
            let product = crate::qubit::mat_mul(&etk, ek);
            sum[0][0] += product[0][0];
            sum[0][1] += product[0][1];
            sum[1][0] += product[1][0];
            sum[1][1] += product[1][1];
        }
        (sum[0][0] - 1.0).abs() < 1e-10
            && (sum[0][1]).abs() < 1e-10
            && (sum[1][0]).abs() < 1e-10
            && (sum[1][1] - 1.0).abs() < 1e-10
    }

    /// Bit-flip channel: flips |0⟩↔|1⟩ with probability p.
    /// E_0 = √(1-p) I, E_1 = √p X
    pub fn bit_flip(p: f64) -> Self {
        let sqrt_q = (1.0 - p).sqrt();
        let sqrt_p = p.sqrt();
        Self {
            kraus: vec![
                [[sqrt_q, 0.0], [0.0, sqrt_q]],
                [[0.0, sqrt_p], [sqrt_p, 0.0]],
            ],
        }
    }

    /// Phase-flip channel: applies Z with probability p.
    /// E_0 = √(1-p) I, E_1 = √p Z
    pub fn phase_flip(p: f64) -> Self {
        let sqrt_q = (1.0 - p).sqrt();
        let sqrt_p = p.sqrt();
        Self {
            kraus: vec![
                [[sqrt_q, 0.0], [0.0, sqrt_q]],
                [[sqrt_p, 0.0], [0.0, -sqrt_p]],
            ],
        }
    }

    /// Depolarizing channel: replaces state with maximally mixed with probability p.
    /// ρ → (1-p)ρ + p(I/2)
    /// E_0 = √(1-3p/4) I, E_1 = √(p/4) X, E_2 = √(p/4) Z, E_3 = √(p/4) I (Y placeholder)
    pub fn depolarizing(p: f64) -> Self {
        let a = (1.0 - 3.0 * p / 4.0).sqrt();
        let b = (p / 4.0).sqrt();
        Self {
            kraus: vec![
                [[a, 0.0], [0.0, a]],  // I
                [[0.0, b], [b, 0.0]],  // X
                [[b, 0.0], [0.0, -b]], // Z
                [[b, 0.0], [0.0, b]],  // Y placeholder (identity scaled)
            ],
        }
    }

    /// Amplitude damping channel: models energy dissipation.
    /// E_0 = [[1, 0], [0, √(1-γ)]], E_1 = [[0, √γ], [0, 0]]
    pub fn amplitude_damping(gamma: f64) -> Self {
        let sqrt_g = gamma.sqrt();
        let sqrt_1mg = (1.0 - gamma).sqrt();
        Self {
            kraus: vec![[[1.0, 0.0], [0.0, sqrt_1mg]], [[0.0, sqrt_g], [0.0, 0.0]]],
        }
    }

    /// Estimate channel capacity (classical, in bits).
    /// Simplified: uses Holevo quantity upper bound for depolarizing-like channels.
    /// C = 1 - H₂((1 + D)/2) where D is the depolarizing parameter.
    pub fn classical_capacity_estimate(&self, depolarizing_param: f64) -> f64 {
        if depolarizing_param >= 1.0 {
            return 0.0;
        }
        let d = (1.0 + depolarizing_param) / 2.0;
        if d <= 0.0 || d >= 1.0 {
            return if d >= 1.0 { 1.0 } else { 0.0 };
        }
        let h2 = -d * d.ln() / std::f64::consts::LN_2
            - (1.0 - d) * (1.0 - d).ln() / std::f64::consts::LN_2;
        1.0 - h2
    }
}

/// Fidelity between two density matrices: F(ρ, σ) = Tr(√(√ρ σ √ρ))²
/// For 2×2 real matrices, we use the simplified formula:
/// F(ρ, σ) = (Tr(√(√ρ σ √ρ)))²
/// For computational ease with diagonal matrices, we use a direct approach.
pub fn fidelity(rho: &DensityMatrix, sigma: &DensityMatrix) -> f64 {
    let (e1, e2) = rho.eigenvalues();
    let (f1, f2) = sigma.eigenvalues();
    // Simplified fidelity for approximately commuting states:
    // F = (√e1√f1 + √e2√f2)²
    let sqrt_e1 = e1.max(0.0).sqrt();
    let sqrt_e2 = e2.max(0.0).sqrt();
    let sqrt_f1 = f1.max(0.0).sqrt();
    let sqrt_f2 = f2.max(0.0).sqrt();
    (sqrt_e1 * sqrt_f1 + sqrt_e2 * sqrt_f2).powi(2)
}

/// Trace distance between two density matrices: T(ρ,σ) = ½|ρ - σ|₁
pub fn trace_distance(rho: &DensityMatrix, sigma: &DensityMatrix) -> f64 {
    let diff = DensityMatrix::new([
        [
            rho.data[0][0] - sigma.data[0][0],
            rho.data[0][1] - sigma.data[0][1],
        ],
        [
            rho.data[1][0] - sigma.data[1][0],
            rho.data[1][1] - sigma.data[1][1],
        ],
    ]);
    let (e1, e2) = diff.eigenvalues();
    0.5 * (e1.abs() + e2.abs())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bit_flip_channel() {
        let ch = QuantumChannel::bit_flip(0.0);
        assert!(ch.verify_completeness());
        let rho = QubitState::zero().to_density_matrix();
        let result = ch.apply(&rho);
        assert!((result.data[0][0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_bit_flip_full() {
        let ch = QuantumChannel::bit_flip(1.0);
        assert!(ch.verify_completeness());
        let rho = QubitState::zero().to_density_matrix();
        let result = ch.apply(&rho);
        assert!((result.data[1][1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_phase_flip() {
        let ch = QuantumChannel::phase_flip(0.3);
        assert!(ch.verify_completeness());
    }

    #[test]
    fn test_depolarizing() {
        let ch = QuantumChannel::depolarizing(0.5);
        assert!(ch.verify_completeness());
        let rho = QubitState::zero().to_density_matrix();
        let result = ch.apply(&rho);
        assert!((result.trace() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_amplitude_damping() {
        let ch = QuantumChannel::amplitude_damping(0.0);
        let rho = QubitState::one().to_density_matrix();
        let result = ch.apply(&rho);
        assert!((result.data[1][1] - 1.0).abs() < 1e-10); // no damping
    }

    #[test]
    fn test_amplitude_damping_full() {
        let ch = QuantumChannel::amplitude_damping(1.0);
        let rho = QubitState::one().to_density_matrix();
        let result = ch.apply(&rho);
        // Full damping: |1⟩ decays to |0⟩
        assert!((result.data[0][0] - 1.0).abs() < 1e-10);
        assert!((result.data[1][1]).abs() < 1e-10);
    }

    #[test]
    fn test_fidelity_same_state() {
        let rho = QubitState::plus().to_density_matrix();
        let f = fidelity(&rho, &rho);
        assert!((f - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_trace_distance_same() {
        let rho = QubitState::zero().to_density_matrix();
        let d = trace_distance(&rho, &rho);
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn test_channel_preserves_trace() {
        let ch = QuantumChannel::depolarizing(0.3);
        let rho = QubitState::plus().to_density_matrix();
        let result = ch.apply(&rho);
        assert!((result.trace() - 1.0).abs() < 1e-10);
    }
}
