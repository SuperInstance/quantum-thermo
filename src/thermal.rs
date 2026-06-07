//! Thermal states, partition functions, and thermodynamic quantities.
//!
//! A thermal (Gibbs) state is ρ_β = e^{-βH} / Z where β = 1/(k_B T)
//! and Z = Tr(e^{-βH}) is the partition function.

use crate::entropy_quantum::von_neumann_entropy;
use crate::qubit::DensityMatrix;
use serde::{Deserialize, Serialize};

/// Inverse temperature β = 1/(k_B T). We use natural units (k_B = 1).
pub type Beta = f64;

/// A thermal state for a 2-level quantum system.
///
/// Hamiltonian H = diag(0, ε) with energy gap ε.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThermalState {
    /// Energy gap ε between ground and excited state
    pub energy_gap: f64,
    /// Inverse temperature β = 1/T (k_B = 1)
    pub beta: f64,
}

impl ThermalState {
    /// Create a new thermal state with given energy gap and inverse temperature.
    pub fn new(energy_gap: f64, beta: f64) -> Self {
        Self { energy_gap, beta }
    }

    /// Create from temperature T (β = 1/T). Returns infinite-β state for T=0.
    pub fn from_temperature(energy_gap: f64, temperature: f64) -> Self {
        let beta = if temperature < 1e-15 {
            1e15 // very cold
        } else {
            1.0 / temperature
        };
        Self { energy_gap, beta }
    }

    /// Partition function Z = Tr(e^{-βH}) = 1 + e^{-βε}
    pub fn partition_function(&self) -> f64 {
        1.0 + (-self.beta * self.energy_gap).exp()
    }

    /// Gibbs state as a density matrix: ρ_β = e^{-βH} / Z
    pub fn density_matrix(&self) -> DensityMatrix {
        let z = self.partition_function();
        let p1 = (-self.beta * self.energy_gap).exp() / z;
        let p0 = 1.0 / z;
        DensityMatrix::new([[p0, 0.0], [0.0, p1]])
    }

    /// Population of ground state: p₀ = 1/Z
    pub fn ground_population(&self) -> f64 {
        1.0 / self.partition_function()
    }

    /// Population of excited state: p₁ = e^{-βε}/Z
    pub fn excited_population(&self) -> f64 {
        (-self.beta * self.energy_gap).exp() / self.partition_function()
    }

    /// Thermal expectation value ⟨H⟩ = ε · p₁
    pub fn mean_energy(&self) -> f64 {
        self.energy_gap * self.excited_population()
    }

    /// Heat capacity C = d⟨H⟩/dT = (βε)² e^{-βε} / (1 + e^{-βε})²
    pub fn heat_capacity(&self) -> f64 {
        let x = self.beta * self.energy_gap;
        let exp_neg = (-x).exp();
        x * x * exp_neg / ((1.0 + exp_neg) * (1.0 + exp_neg))
    }

    /// Thermal entropy S = -Σ pᵢ ln(pᵢ)
    pub fn entropy(&self) -> f64 {
        let rho = self.density_matrix();
        von_neumann_entropy(&rho)
    }

    /// Free energy F = -T ln(Z) = -(1/β) ln(Z)
    pub fn free_energy(&self) -> f64 {
        -self.partition_function().ln() / self.beta
    }

    /// Internal energy U = ⟨H⟩
    pub fn internal_energy(&self) -> f64 {
        self.mean_energy()
    }
}

/// Compute the partition function for a 2-level system at given β and ε.
pub fn partition_function_2level(beta: f64, energy_gap: f64) -> f64 {
    1.0 + (-beta * energy_gap).exp()
}

/// Boltzmann factor e^{-β E}
pub fn boltzmann_factor(beta: f64, energy: f64) -> f64 {
    (-beta * energy).exp()
}

/// Mean energy of a 2-level system: ⟨E⟩ = ε / (1 + e^{βε})
pub fn mean_energy_2level(beta: f64, energy_gap: f64) -> f64 {
    energy_gap / (1.0 + (beta * energy_gap).exp())
}

/// Relative population of excited state.
pub fn boltzmann_population(beta: f64, energy_gap: f64) -> f64 {
    1.0 / (1.0 + (beta * energy_gap).exp())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_partition_function() {
        let ts = ThermalState::new(1.0, 1.0);
        let z = ts.partition_function();
        let expected = 1.0 + (-1.0_f64).exp();
        assert!((z - expected).abs() < 1e-10);
    }

    #[test]
    fn test_ground_state_zero_temp() {
        let ts = ThermalState::new(1.0, 100.0); // very cold
        assert!((ts.ground_population() - 1.0).abs() < 1e-6);
        assert!(ts.excited_population() < 1e-10);
    }

    #[test]
    fn test_equal_population_high_temp() {
        let ts = ThermalState::new(1.0, 0.001); // very hot
        assert!((ts.ground_population() - 0.5).abs() < 1e-3);
        assert!((ts.excited_population() - 0.5).abs() < 1e-3);
    }

    #[test]
    fn test_mean_energy() {
        let ts = ThermalState::new(1.0, 1.0);
        let e = ts.mean_energy();
        assert!(e > 0.0 && e < 1.0);
    }

    #[test]
    fn test_heat_capacity() {
        let ts = ThermalState::new(1.0, 1.0);
        let c = ts.heat_capacity();
        assert!(c > 0.0);
    }

    #[test]
    fn test_free_energy() {
        let ts = ThermalState::new(1.0, 1.0);
        let f = ts.free_energy();
        assert!(f < 0.0); // F = -ln(Z)/β, Z > 1, so F < 0
    }

    #[test]
    fn test_from_temperature() {
        let ts = ThermalState::from_temperature(2.0, 1.0);
        assert!((ts.beta - 1.0).abs() < 1e-10);
        assert!((ts.energy_gap - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_density_matrix_trace() {
        let ts = ThermalState::new(1.0, 0.5);
        let rho = ts.density_matrix();
        assert!((rho.trace() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_entropy_high_temp() {
        let ts = ThermalState::new(1.0, 0.001);
        let s = ts.entropy();
        // At high T, should approach ln(2)
        assert!((s - std::f64::consts::LN_2).abs() < 0.01);
    }
}
