//! Quantum work extraction via two-point measurement scheme.
//!
//! Work is defined as W = E_final - E_initial where energies are obtained
//! from projective measurements before and after the process.

use serde::{Deserialize, Serialize};

/// A single work trajectory from the two-point measurement scheme.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WorkTrajectory {
    /// Initial energy eigenvalue
    pub e_initial: f64,
    /// Final energy eigenvalue
    pub e_final: f64,
    /// Probability of this trajectory
    pub probability: f64,
}

impl WorkTrajectory {
    /// Work done in this trajectory: W = E_final - E_initial
    pub fn work(&self) -> f64 {
        self.e_final - self.e_initial
    }
}

/// Result of a work distribution analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WorkDistribution {
    /// Individual work trajectories
    pub trajectories: Vec<WorkTrajectory>,
}

impl WorkDistribution {
    /// Create from a list of trajectories.
    pub fn new(trajectories: Vec<WorkTrajectory>) -> Self {
        Self { trajectories }
    }

    /// Compute the two-point work distribution for a 2-level system.
    ///
    /// Given initial thermal state at β with Hamiltonian H_i = diag(0, ε_i)
    /// and final Hamiltonian H_f = diag(0, ε_f), compute all trajectories.
    pub fn two_level(beta: f64, energy_initial: f64, energy_final: f64) -> Self {
        let z = 1.0 + (-beta * energy_initial).exp();
        let p0 = 1.0 / z;
        let p1 = (-beta * energy_initial).exp() / z;

        let trajectories = vec![
            WorkTrajectory {
                e_initial: 0.0,
                e_final: 0.0,
                probability: p0,
            },
            WorkTrajectory {
                e_initial: energy_initial,
                e_final: energy_final,
                probability: p1,
            },
        ];

        Self { trajectories }
    }

    /// Compute from explicit transition probabilities.
    /// trans_probs[i][j] = probability of measuring E_f=j | E_i=i
    pub fn from_transitions_2d(
        beta: f64,
        energies_initial: &[f64],
        energies_final: &[f64],
        trans_probs: &[[f64; 2]; 2],
    ) -> Self {
        let z = energies_initial
            .iter()
            .map(|e| (-beta * e).exp())
            .sum::<f64>();
        let mut trajectories = Vec::new();
        for (i, &ei) in energies_initial.iter().enumerate() {
            let pi = (-beta * ei).exp() / z;
            for (j, &ef) in energies_final.iter().enumerate() {
                let pj = trans_probs[i][j];
                let prob = pi * pj;
                if prob > 1e-15 {
                    trajectories.push(WorkTrajectory {
                        e_initial: ei,
                        e_final: ef,
                        probability: prob,
                    });
                }
            }
        }
        Self { trajectories }
    }

    /// Mean work ⟨W⟩
    pub fn mean_work(&self) -> f64 {
        self.trajectories
            .iter()
            .map(|t| t.probability * t.work())
            .sum()
    }

    /// Work variance Var(W) = ⟨W²⟩ - ⟨W⟩²
    pub fn work_variance(&self) -> f64 {
        let mean = self.mean_work();
        let mean2: f64 = self
            .trajectories
            .iter()
            .map(|t| t.probability * t.work().powi(2))
            .sum();
        mean2 - mean * mean
    }

    /// Jarzynski equality: ⟨e^{-βW}⟩ should equal Z_f / Z_i = 1 for equilibrium.
    pub fn jarzynski_average(&self, beta: f64) -> f64 {
        self.trajectories
            .iter()
            .map(|t| t.probability * (-beta * t.work()).exp())
            .sum()
    }

    /// Verify the Jarzynski equality: ⟨e^{-βW}⟩ = Z_f / Z_i
    pub fn verify_jarzynski(&self, beta: f64, z_initial: f64, z_final: f64) -> bool {
        let avg = self.jarzynski_average(beta);
        let expected = z_final / z_initial;
        (avg - expected).abs() < 0.01 // allow numerical tolerance
    }
}

/// Fluctuation theorem: P(W)/P(-W) = e^{βW}
/// Returns the ratio P(W)/P(-W) for given work value.
pub fn fluctuation_ratio(work: f64, beta: f64) -> f64 {
    (beta * work).exp()
}

/// Compute Jarzynski equality for a parameter change in 2-level system.
/// Returns ⟨e^{-βW}⟩ = Z_f(β) / Z_i(β)
pub fn jarzynski_equality(beta: f64, energy_initial: f64, energy_final: f64) -> f64 {
    let zi = 1.0 + (-beta * energy_initial).exp();
    let zf = 1.0 + (-beta * energy_final).exp();
    zf / zi
}

/// Crooks fluctuation theorem: P_fwd(W) / P_rev(-W) = e^{β(W - ΔF)}
pub fn crooks_ratio(work: f64, beta: f64, delta_f: f64) -> f64 {
    (beta * (work - delta_f)).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_work_trajectory() {
        let t = WorkTrajectory {
            e_initial: 0.0,
            e_final: 1.0,
            probability: 0.5,
        };
        assert!((t.work() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_two_level_distribution() {
        let dist = WorkDistribution::two_level(1.0, 1.0, 2.0);
        let total_prob: f64 = dist.trajectories.iter().map(|t| t.probability).sum();
        assert!((total_prob - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_mean_work() {
        let dist = WorkDistribution::two_level(1.0, 1.0, 2.0);
        let w = dist.mean_work();
        // W = (0→0)·p0 + (1→2)·p1 = 0·p0 + 1·p1
        let p1 = 1.0 / (1.0 + 1.0_f64.exp());
        assert!((w - p1).abs() < 1e-10);
    }

    #[test]
    fn test_jarzynski_average() {
        let dist = WorkDistribution::two_level(1.0, 1.0, 2.0);
        let avg = dist.jarzynski_average(1.0);
        // Should equal Z_f/Z_i
        let zi = 1.0 + (-1.0_f64).exp();
        let zf = 1.0 + (-2.0_f64).exp();
        assert!((avg - zf / zi).abs() < 1e-10);
    }

    #[test]
    fn test_jarzynski_identity() {
        // When no protocol (same initial and final), ⟨e^{-βW}⟩ = 1
        let dist = WorkDistribution::two_level(1.0, 1.0, 1.0);
        let avg = dist.jarzynski_average(1.0);
        assert!((avg - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_fluctuation_ratio() {
        let r = fluctuation_ratio(1.0, 1.0);
        assert!((r - std::f64::consts::E).abs() < 1e-10);
    }

    #[test]
    fn test_jarzynski_equality_fn() {
        let val = jarzynski_equality(1.0, 1.0, 2.0);
        let zi = 1.0 + (-1.0_f64).exp();
        let zf = 1.0 + (-2.0_f64).exp();
        assert!((val - zf / zi).abs() < 1e-10);
    }

    #[test]
    fn test_work_variance() {
        let dist = WorkDistribution::two_level(1.0, 1.0, 2.0);
        let v = dist.work_variance();
        assert!(v >= 0.0);
    }

    #[test]
    fn test_from_transitions() {
        let trans: [[f64; 2]; 2] = [
            [1.0, 0.0], // initial level 0 → stays 0
            [0.0, 1.0], // initial level 1 → stays 1
        ];
        let dist = WorkDistribution::from_transitions_2d(1.0, &[0.0, 1.0], &[0.0, 2.0], &trans);
        assert_eq!(dist.trajectories.len(), 2);
    }
}
