//! # quantum-thermo
//!
//! Quantum thermodynamics: quantum limits on computation, heat, and information.
//!
//! This crate provides tools for simulating quantum thermodynamic systems,
//! including qubit states, quantum channels, thermal states, work extraction,
//! von Neumann entropy, and resource theories.
//!
//! ## Modules
//!
//! - [`qubit`] — Qubit states, Pauli matrices, density matrices, Bloch sphere
//! - [`channel`] — Quantum channels (Kraus representation)
//! - [`thermal`] — Thermal/Gibbs states, partition functions
//! - [`work`] — Quantum work (two-point measurement scheme)
//! - [`entropy_quantum`] — Von Neumann entropy and related quantities
//! - [`resource`] — Quantum resource theories

pub mod channel;
pub mod entropy_quantum;
pub mod qubit;
pub mod resource;
pub mod thermal;
pub mod work;
