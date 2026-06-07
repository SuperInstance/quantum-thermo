//! Quantum resource theories: entanglement, athermality, and coherence.
//!
//! A resource theory defines (1) a set of free states and (2) a set of free
//! operations. States outside the free set are resources. This module provides
//! tools for quantifying and manipulating quantum resources.

use crate::entropy_quantum::{BipartiteState, entanglement_entropy, von_neumann_entropy};
use crate::qubit::DensityMatrix;
use serde::{Deserialize, Serialize};

/// The type of resource theory.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum ResourceTheoryType {
    /// Entanglement theory: free states = separable, free ops = LOCC
    Entanglement,
    /// Thermal resource theory: free states = Gibbs states, free ops = thermal ops
    Athermality,
    /// Coherence theory: free states = diagonal, free ops = incoherent
    Coherence,
}

/// A resource monotone value.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResourceValue {
    /// Name of the monotone
    pub name: String,
    /// Value
    pub value: f64,
}

impl ResourceValue {
    pub fn new(name: &str, value: f64) -> Self {
        Self {
            name: name.to_string(),
            value,
        }
    }
}

/// Check if a single-qubit state is a valid free state for the given resource theory.
pub fn is_free_state(rho: &DensityMatrix, theory: ResourceTheoryType) -> bool {
    match theory {
        ResourceTheoryType::Entanglement => {
            // Single-qubit states are always separable (no entanglement possible)
            true
        }
        ResourceTheoryType::Athermality => {
            // Thermal states are free. For simplicity, check if diagonal.
            rho.data[0][1].abs() < 1e-10 && rho.data[1][0].abs() < 1e-10
        }
        ResourceTheoryType::Coherence => {
            // Free states are diagonal in the computational basis
            rho.data[0][1].abs() < 1e-10 && rho.data[1][0].abs() < 1e-10
        }
    }
}

/// Check if a bipartite state is separable (free in entanglement theory).
/// A state is separable if S(ρ_A) ≥ entanglement_entropy.
pub fn is_separable(state_ab: &BipartiteState) -> bool {
    let rho_a = state_ab.partial_trace_b();
    let _s_a = von_neumann_entropy(&rho_a);
    let ee = entanglement_entropy(state_ab);
    // For pure states: separable iff entanglement entropy = 0
    ee.abs() < 1e-10
}

/// Compute resource monotones for a single-qubit state.
pub fn resource_monotones(rho: &DensityMatrix, theory: ResourceTheoryType) -> Vec<ResourceValue> {
    match theory {
        ResourceTheoryType::Entanglement => {
            // Single qubit can't have entanglement
            vec![ResourceValue::new("entanglement", 0.0)]
        }
        ResourceTheoryType::Athermality => {
            // Relative entropy of athermality: S(ρ || ρ_thermal)
            let thermal = DensityMatrix::maximally_mixed();
            let s_rho = von_neumann_entropy(rho);
            let s_thermal = von_neumann_entropy(&thermal);
            vec![
                ResourceValue::new("athermality", s_thermal - s_rho),
                ResourceValue::new("purity", rho.purity()),
            ]
        }
        ResourceTheoryType::Coherence => {
            // L1-norm of coherence: sum of absolute values of off-diagonal elements
            let l1 = rho.data[0][1].abs() + rho.data[1][0].abs();
            // Relative entropy of coherence: S(ρ_diag) - S(ρ)
            let diag = DensityMatrix::new([[rho.data[0][0], 0.0], [0.0, rho.data[1][1]]]);
            let coh = von_neumann_entropy(&diag) - von_neumann_entropy(rho);
            vec![
                ResourceValue::new("l1_coherence", l1),
                ResourceValue::new("relative_entropy_coherence", coh),
            ]
        }
    }
}

/// Entanglement monotones for a bipartite state.
pub fn entanglement_monotones(state_ab: &BipartiteState) -> Vec<ResourceValue> {
    let ee = entanglement_entropy(state_ab);
    let mi = crate::entropy_quantum::mutual_information(state_ab);
    vec![
        ResourceValue::new("entanglement_entropy", ee),
        ResourceValue::new("mutual_information", mi),
    ]
}

/// Result of a state conversion check.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConversionResult {
    /// Whether the conversion is possible
    pub possible: bool,
    /// Reason
    pub reason: String,
}

/// Check if state ρ can be converted to σ using free operations.
pub fn can_convert(
    rho: &DensityMatrix,
    sigma: &DensityMatrix,
    theory: ResourceTheoryType,
) -> ConversionResult {
    match theory {
        ResourceTheoryType::Entanglement => {
            // Nielsen's theorem: ρ → σ by LOCC iff λ_ρ ≻ λ_σ (majorization)
            ConversionResult {
                possible: true,
                reason: "Single-qubit states are always interconvertible by LOCC (no entanglement)"
                    .to_string(),
            }
        }
        ResourceTheoryType::Athermality => {
            // Check if resource monotone decreases
            let mono_rho = resource_monotones(rho, theory);
            let mono_sigma = resource_monotones(sigma, theory);
            if mono_rho[0].value >= mono_sigma[0].value {
                ConversionResult {
                    possible: true,
                    reason: format!(
                        "Athermality decreases: {:.4} → {:.4}",
                        mono_rho[0].value, mono_sigma[0].value
                    ),
                }
            } else {
                ConversionResult {
                    possible: false,
                    reason: format!(
                        "Athermality would increase: {:.4} → {:.4}",
                        mono_rho[0].value, mono_sigma[0].value
                    ),
                }
            }
        }
        ResourceTheoryType::Coherence => {
            let coh_rho = resource_monotones(rho, theory);
            let coh_sigma = resource_monotones(sigma, theory);
            if coh_rho[0].value >= coh_sigma[0].value {
                ConversionResult {
                    possible: true,
                    reason: format!(
                        "Coherence decreases: {:.4} → {:.4}",
                        coh_rho[0].value, coh_sigma[0].value
                    ),
                }
            } else {
                ConversionResult {
                    possible: false,
                    reason: format!(
                        "Coherence would increase: {:.4} → {:.4}",
                        coh_rho[0].value, coh_sigma[0].value
                    ),
                }
            }
        }
    }
}

/// Simulate a simple LOCC operation: local measurement on subsystem A.
/// Returns the post-measurement bipartite state.
pub fn locc_measure_a(state_ab: &BipartiteState, outcome: usize) -> BipartiteState {
    assert!(outcome < 2);
    let d = &state_ab.data;
    let idx = |r: usize, c: usize| r * 4 + c;

    // Measure A in computational basis
    // If outcome = 0: keep rows/cols 0,1 (A in |0⟩)
    // If outcome = 1: keep rows/cols 2,3 (A in |1⟩)
    let mut result = vec![0.0; 16];
    let offset = outcome * 2;
    let mut norm = 0.0;

    for i in 0..2 {
        for j in 0..2 {
            let val = d[idx(offset + i, offset + j)];
            result[idx(i, j)] = val;
            norm += val;
        }
    }

    // Normalize
    if norm > 1e-15 {
        for val in result.iter_mut() {
            *val /= norm;
        }
    }

    BipartiteState::new(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qubit::QubitState;

    #[test]
    fn test_free_state_entanglement() {
        let rho = QubitState::zero().to_density_matrix();
        assert!(is_free_state(&rho, ResourceTheoryType::Entanglement));
    }

    #[test]
    fn test_free_state_coherence() {
        let rho = QubitState::zero().to_density_matrix();
        assert!(is_free_state(&rho, ResourceTheoryType::Coherence));
    }

    #[test]
    fn test_not_free_coherence() {
        let rho = QubitState::plus().to_density_matrix();
        assert!(!is_free_state(&rho, ResourceTheoryType::Coherence));
    }

    #[test]
    fn test_coherence_monotone() {
        let rho = QubitState::plus().to_density_matrix();
        let monos = resource_monotones(&rho, ResourceTheoryType::Coherence);
        assert!(monos[0].value > 0.0); // L1 coherence > 0
    }

    #[test]
    fn test_bell_entanglement() {
        let state = BipartiteState::bell_state();
        let monos = entanglement_monotones(&state);
        assert!(monos[0].value > 0.0);
    }

    #[test]
    fn test_product_separable() {
        let state = BipartiteState::product_00();
        assert!(is_separable(&state));
    }

    #[test]
    fn test_bell_not_separable() {
        let state = BipartiteState::bell_state();
        assert!(!is_separable(&state));
    }

    #[test]
    fn test_can_convert_coherence() {
        let rho = QubitState::plus().to_density_matrix();
        let sigma = QubitState::zero().to_density_matrix();
        let result = can_convert(&rho, &sigma, ResourceTheoryType::Coherence);
        assert!(result.possible);
    }

    #[test]
    fn test_cannot_increase_coherence() {
        let rho = QubitState::zero().to_density_matrix();
        let sigma = QubitState::plus().to_density_matrix();
        let result = can_convert(&rho, &sigma, ResourceTheoryType::Coherence);
        assert!(!result.possible);
    }

    #[test]
    fn test_locc_measure() {
        let state = BipartiteState::bell_state();
        let post = locc_measure_a(&state, 0);
        // After measuring A=0, B should be in |0⟩
        assert!((post.data[0] - 1.0).abs() < 1e-10);
    }
}
