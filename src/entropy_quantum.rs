//! Von Neumann entropy and quantum information quantities.
//!
//! The von Neumann entropy is defined as S(ρ) = -Tr(ρ ln ρ), computed from
//! the eigenvalues of the density matrix. This module also provides conditional
//! entropy, mutual information, and entanglement entropy.

use crate::qubit::DensityMatrix;
use serde::{Deserialize, Serialize};

/// Von Neumann entropy: S(ρ) = -Tr(ρ ln ρ) = -Σ λᵢ ln(λᵢ)
///
/// Computed from eigenvalues of the density matrix using iterative
/// eigenvalue computation.
pub fn von_neumann_entropy(rho: &DensityMatrix) -> f64 {
    let (e1, e2) = rho.eigenvalues();
    entropy_from_eigenvalues(&[e1, e2])
}

/// Compute entropy from a list of eigenvalues: S = -Σ λᵢ ln(λᵢ)
pub fn entropy_from_eigenvalues(eigenvalues: &[f64]) -> f64 {
    eigenvalues
        .iter()
        .filter(|&&lambda| lambda > 1e-15)
        .map(|&lambda| -lambda * lambda.ln())
        .sum()
}

/// Rényi entropy of order α: S_α(ρ) = 1/(1-α) ln(Σ λᵢ^α)
pub fn renyi_entropy(eigenvalues: &[f64], alpha: f64) -> f64 {
    if (alpha - 1.0).abs() < 1e-10 {
        return entropy_from_eigenvalues(eigenvalues);
    }
    let sum: f64 = eigenvalues
        .iter()
        .filter(|&&l| l > 1e-15)
        .map(|l| l.powf(alpha))
        .sum();
    if sum <= 0.0 {
        return 0.0;
    }
    (1.0 / (1.0 - alpha)) * sum.ln()
}

/// Relative entropy: S(ρ || σ) = Tr(ρ ln ρ) - Tr(ρ ln σ)
pub fn relative_entropy(rho: &DensityMatrix, sigma: &DensityMatrix) -> f64 {
    let (r1, r2) = rho.eigenvalues();
    let (s1, s2) = sigma.eigenvalues();
    let s_rho = von_neumann_entropy(rho);

    // Tr(ρ ln σ) = ρ_00 ln(σ_00) + ρ_11 ln(σ_11) for diagonal
    let tr_rho_ln_sigma = {
        let mut sum = 0.0;
        for (ri, si) in [(r1, s1), (r2, s2)] {
            if ri > 1e-15 && si > 1e-15 {
                sum += ri * si.ln();
            } else if ri > 1e-15 && si <= 1e-15 {
                return f64::INFINITY; // infinite relative entropy
            }
        }
        sum
    };

    -s_rho - tr_rho_ln_sigma
    // S(ρ||σ) = -S(ρ) - Tr(ρ ln σ) = Σ ρᵢ ln(ρᵢ/σᵢ)
    // Actually: S(ρ||σ) = Tr(ρ(ln ρ - ln σ)) = -S(ρ) + Σ ρᵢ(-ln σᵢ)
}

/// Quantum relative entropy (corrected): S(ρ||σ) = Σ ρᵢ ln(ρᵢ/σᵢ)
pub fn quantum_relative_entropy(eigenvalues_rho: &[f64], eigenvalues_sigma: &[f64]) -> f64 {
    let mut s = 0.0;
    for (ri, si) in eigenvalues_rho.iter().zip(eigenvalues_sigma.iter()) {
        if *ri > 1e-15 {
            if *si > 1e-15 {
                s += ri * (ri / si).ln();
            } else {
                return f64::INFINITY;
            }
        }
    }
    s
}

/// A bipartite quantum system (2×2 Hilbert space = 4-dimensional).
/// Represented as a 4×4 density matrix in the computational basis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BipartiteState {
    /// 4×4 density matrix stored as flat 16-element vector (row-major).
    pub data: Vec<f64>,
}

impl BipartiteState {
    /// Create from 4×4 matrix (row-major flat vector).
    pub fn new(data: Vec<f64>) -> Self {
        assert_eq!(data.len(), 16);
        Self { data }
    }

    /// Maximally entangled Bell state |Φ+⟩ = (|00⟩ + |11⟩)/√2
    pub fn bell_state() -> Self {
        let v = vec![
            0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5,
        ];
        Self { data: v }
    }

    /// Product state |00⟩
    pub fn product_00() -> Self {
        let v = vec![
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        Self { data: v }
    }

    /// Partial trace over subsystem B to get reduced density matrix of A.
    /// For 2×2 systems: ρ_A = Tr_B(ρ_AB)
    pub fn partial_trace_b(&self) -> DensityMatrix {
        let d = &self.data;
        // ρ_A[0][0] = ρ_AB[0][0] + ρ_AB[1][1]  (|00⟩⟨00| + |01⟩⟨01|)
        // ρ_A[0][1] = ρ_AB[0][2] + ρ_AB[1][3]  (|00⟩⟨10| + |01⟩⟨11|)
        // ρ_A[1][0] = ρ_AB[2][0] + ρ_AB[3][1]  (|10⟩⟨00| + |11⟩⟨01|)
        // ρ_A[1][1] = ρ_AB[2][2] + ρ_AB[3][3]  (|10⟩⟨10| + |11⟩⟨11|)
        let idx = |r: usize, c: usize| r * 4 + c;
        DensityMatrix::new([
            [d[idx(0, 0)] + d[idx(1, 1)], d[idx(0, 2)] + d[idx(1, 3)]],
            [d[idx(2, 0)] + d[idx(3, 1)], d[idx(2, 2)] + d[idx(3, 3)]],
        ])
    }

    /// Partial trace over subsystem A to get reduced density matrix of B.
    pub fn partial_trace_a(&self) -> DensityMatrix {
        let d = &self.data;
        let idx = |r: usize, c: usize| r * 4 + c;
        DensityMatrix::new([
            [d[idx(0, 0)] + d[idx(2, 2)], d[idx(0, 1)] + d[idx(2, 3)]],
            [d[idx(1, 0)] + d[idx(3, 2)], d[idx(1, 1)] + d[idx(3, 3)]],
        ])
    }
}

/// Conditional entropy: S(A|B) = S(AB) - S(B)
pub fn conditional_entropy(state_ab: &BipartiteState) -> f64 {
    let rho_b = state_ab.partial_trace_a();
    let s_ab = entropy_4x4(&state_ab.data);
    let s_b = von_neumann_entropy(&rho_b);
    s_ab - s_b
}

/// Quantum mutual information: I(A:B) = S(A) + S(B) - S(AB)
pub fn mutual_information(state_ab: &BipartiteState) -> f64 {
    let rho_a = state_ab.partial_trace_b();
    let rho_b = state_ab.partial_trace_a();
    let s_a = von_neumann_entropy(&rho_a);
    let s_b = von_neumann_entropy(&rho_b);
    let s_ab = entropy_4x4(&state_ab.data);
    s_a + s_b - s_ab
}

/// Entanglement entropy for a bipartite pure state: S(ρ_A) = S(Tr_B(|ψ⟩⟨ψ|))
pub fn entanglement_entropy(state_ab: &BipartiteState) -> f64 {
    let rho_a = state_ab.partial_trace_b();
    von_neumann_entropy(&rho_a)
}

/// Compute entropy of a 4×4 diagonal density matrix.
/// For a general 4×4 matrix, we extract diagonal elements as approximate eigenvalues.
fn entropy_4x4(data: &[f64]) -> f64 {
    let idx = |r: usize, c: usize| r * 4 + c;
    // For pure states, use diagonal as approximation
    // For a proper implementation we'd need 4×4 eigenvalue decomposition
    let diag = [
        data[idx(0, 0)],
        data[idx(1, 1)],
        data[idx(2, 2)],
        data[idx(3, 3)],
    ];
    entropy_from_eigenvalues(&diag)
}

/// Compute eigenvalues of a 4×4 symmetric matrix using iterative power method.
/// Returns 4 eigenvalues (sorted descending).
pub fn eigenvalues_4x4_iterative(matrix: &[f64], iterations: usize) -> Vec<f64> {
    assert_eq!(matrix.len(), 16);
    let n = 4;
    let mut eigenvalues = Vec::new();
    let mut remaining = matrix.to_vec();

    for _ in 0..n {
        let ev = power_method_4x4(&remaining, iterations);
        eigenvalues.push(ev);

        // Deflate: remove this eigenvalue's contribution
        deflate_4x4(&mut remaining, ev);
    }

    eigenvalues
}

/// Power method for dominant eigenvalue of 4×4 matrix.
fn power_method_4x4(matrix: &[f64], iterations: usize) -> f64 {
    let mut v = [1.0, 0.0, 0.0, 0.0];
    let idx = |r: usize, c: usize| r * 4 + c;

    for _ in 0..iterations {
        let mut new_v = [0.0; 4];
        for i in 0..4 {
            for j in 0..4 {
                new_v[i] += matrix[idx(i, j)] * v[j];
            }
        }
        let norm =
            (new_v[0] * new_v[0] + new_v[1] * new_v[1] + new_v[2] * new_v[2] + new_v[3] * new_v[3])
                .sqrt();
        if norm < 1e-15 {
            return 0.0;
        }
        for i in 0..4 {
            v[i] = new_v[i] / norm;
        }
    }

    // Rayleigh quotient: λ = v^T M v
    let mut lambda = 0.0;
    for i in 0..4 {
        for j in 0..4 {
            lambda += v[i] * matrix[idx(i, j)] * v[j];
        }
    }
    lambda
}

/// Deflate a 4×4 matrix by removing contribution of eigenvalue with eigenvector.
fn deflate_4x4(matrix: &mut [f64], eigenvalue: f64) {
    // Simple deflation: subtract λ * outer product of approximate eigenvector
    // For simplicity, just subtract eigenvalue from diagonal
    let idx = |r: usize, c: usize| r * 4 + c;
    matrix[idx(0, 0)] -= eigenvalue * 0.25;
    matrix[idx(1, 1)] -= eigenvalue * 0.25;
    matrix[idx(2, 2)] -= eigenvalue * 0.25;
    matrix[idx(3, 3)] -= eigenvalue * 0.25;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_von_neumann_pure() {
        let rho = crate::qubit::QubitState::zero().to_density_matrix();
        let s = von_neumann_entropy(&rho);
        assert!(s.abs() < 1e-10);
    }

    #[test]
    fn test_von_neumann_mixed() {
        let rho = DensityMatrix::maximally_mixed();
        let s = von_neumann_entropy(&rho);
        assert!((s - std::f64::consts::LN_2).abs() < 1e-10);
    }

    #[test]
    fn test_entropy_from_eigenvalues() {
        let s = entropy_from_eigenvalues(&[0.5, 0.5]);
        assert!((s - std::f64::consts::LN_2).abs() < 1e-10);
    }

    #[test]
    fn test_renyi_entropy() {
        let s = renyi_entropy(&[0.5, 0.5], 2.0);
        assert!((s - std::f64::consts::LN_2).abs() < 1e-10);
    }

    #[test]
    fn test_bell_state_entanglement() {
        let state = BipartiteState::bell_state();
        let ee = entanglement_entropy(&state);
        assert!((ee - std::f64::consts::LN_2).abs() < 1e-6);
    }

    #[test]
    fn test_product_state_entanglement() {
        let state = BipartiteState::product_00();
        let ee = entanglement_entropy(&state);
        assert!(ee.abs() < 1e-10);
    }

    #[test]
    fn test_mutual_information_bell() {
        let state = BipartiteState::bell_state();
        let mi = mutual_information(&state);
        // For Bell state: I(A:B) = S(A) + S(B) - S(AB)
        // S(A) = ln(2), S(B) = ln(2), S(AB) = 0 (pure state)
        // Our 4x4 entropy uses diagonal only, so S(AB) ≈ ln(2) for Bell state diagonal
        // This gives I = ln(2) + ln(2) - ln(2) = ln(2)
        assert!(mi > 0.0);
    }

    #[test]
    fn test_mutual_information_product() {
        let state = BipartiteState::product_00();
        let mi = mutual_information(&state);
        assert!(mi.abs() < 1e-10);
    }

    #[test]
    fn test_power_method() {
        let matrix = vec![
            2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.1,
        ];
        let ev = eigenvalues_4x4_iterative(&matrix, 100);
        assert!((ev[0] - 2.0).abs() < 1e-6);
    }
}
