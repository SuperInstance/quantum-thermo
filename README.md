# quantum-thermo

**Quantum thermodynamics: quantum limits on computation, heat, and information.**

A pure-Rust library for simulating quantum thermodynamic systems — from qubit states
on the Bloch sphere to fluctuation theorems and entanglement entropy. Built for
physicists, students, and researchers who want concrete, runnable code instead of
abstract promises.

```toml
[dependencies]
quantum-thermo = "0.1"
```

---

## Table of Contents

- [Overview](#overview)
- [Theory](#theory)
  - [Density Matrices](#density-matrices)
  - [Quantum Channels](#quantum-channels-theory)
  - [Thermal States](#thermal-states-theory)
  - [Quantum Work](#quantum-work-theory)
  - [Von Neumann Entropy](#von-neumann-entropy-theory)
  - [Resource Theories](#resource-theories-theory)
- [Architecture](#architecture)
- [Module Reference](#module-reference)
- [Examples](#examples)
  - [1. Bloch Sphere Exploration](#example-1-bloch-sphere-exploration)
  - [2. Quantum Channels and Fidelity](#example-2-quantum-channels-and-fidelity)
  - [3. Thermal States at Various Temperatures](#example-3-thermal-states-at-various-temperatures)
  - [4. Von Neumann Entropy and Entanglement](#example-4-von-neumann-entropy-and-entanglement)
- [Performance Analysis](#performance-analysis)
- [Design Decisions](#design-decisions)
- [Classical vs Quantum Comparison](#classical-vs-quantum-comparison)
- [Glossary](#glossary)
- [References](#references)
- [License](#license)

---

## Overview

Quantum thermodynamics sits at the intersection of quantum mechanics, information theory,
and statistical mechanics. It asks fundamental questions:

- **How much energy does erasing a bit cost?** (Landauer's principle)
- **What are the thermodynamic costs of quantum computation?**
- **Can quantum coherence be used as a thermodynamic resource?**
- **Do fluctuation theorems hold in quantum systems?**

This library provides the computational tools to explore these questions concretely.
Every type is serializable, every function is tested, and every module is documented.

### What This Crate Provides

| Module | What It Does | Key Types |
|--------|-------------|-----------|
| `qubit` | Qubit states, Pauli matrices, density matrices, Bloch sphere | `QubitState`, `DensityMatrix` |
| `channel` | Quantum channels (Kraus representation), noise models | `QuantumChannel` |
| `thermal` | Gibbs states, partition functions, heat capacity | `ThermalState` |
| `work` | Two-point measurement scheme, Jarzynski equality | `WorkDistribution` |
| `entropy_quantum` | Von Neumann entropy, mutual information, entanglement | `BipartiteState` |
| `resource` | Resource theories (entanglement, coherence, athermality) | `ResourceTheoryType` |

---

## Theory

This section provides the mathematical foundations. We use standard physics notation
throughout, with concrete definitions that map directly to the library's types.

### Density Matrices

The density matrix is the fundamental object in quantum statistical mechanics. For an
ensemble of pure states {|ψᵢ⟩} with probabilities {pᵢ}:

```
ρ = Σᵢ pᵢ |ψᵢ⟩⟨ψᵢ|
```

**Properties:**
- **Hermitian**: ρ = ρ†
- **Positive semidefinite**: ⟨φ|ρ|φ⟩ ≥ 0 for all |φ⟩
- **Unit trace**: Tr(ρ) = 1
- **Eigenvalues**: λᵢ ∈ [0, 1], Σ λᵢ = 1

For a single qubit in our real-valued representation, a pure state |ψ⟩ = α|0⟩ + β|1⟩
has density matrix:

```
ρ = |α|²    αβ |
    |βα    |β|²|
```

The **purity** Tr(ρ²) distinguishes pure states (purity = 1) from mixed states (purity < 1).

The **Bloch sphere** parametrizes qubit states as:

```
ρ = ½(I + r⃗·σ⃗)
```

where r⃗ = (rx, ry, rz) is the Bloch vector with |r⃗| ≤ 1, and σ⃗ = (X, Y, Z) are the
Pauli matrices.

**Pauli matrices:**

```
X = |0 1|    Y = |0 -i|    Z = |1  0|
    |1 0|        |i  0|        |0 -1|
```

### Quantum Channels {#quantum-channels-theory}

A quantum channel is a completely positive trace-preserving (CPTP) map. In the **Kraus
representation**, any quantum channel can be written as:

```
ρ' = Σₖ Eₖ ρ Eₖ†
```

where {Eₖ} are the Kraus operators satisfying the completeness relation:

```
Σₖ Eₖ† Eₖ = I
```

**Standard channels implemented:**

| Channel | Kraus Operators | Physical Meaning |
|---------|----------------|------------------|
| BitFlip(p) | √(1-p)I, √p·X | Classical bit flip with probability p |
| PhaseFlip(p) | √(1-p)I, √p·Z | Phase error with probability p |
| Depolarizing(p) | Mixed set | Isotropic noise toward I/2 |
| AmplitudeDamping(γ) | [[1,0],[0,√(1-γ)]], [[0,√γ],[0,0]] | Spontaneous emission |

The **channel capacity** quantifies the maximum rate of classical information transmission.
For a depolarizing channel with parameter p, the Holevo bound gives:

```
C = 1 - H₂((1+p)/2)
```

where H₂ is the binary entropy function.

### Thermal States {#thermal-states-theory}

A system in thermal equilibrium at temperature T is described by the **Gibbs state**:

```
ρ_β = e^{-βH} / Z
```

where β = 1/(k_B T) is the inverse temperature and Z is the **partition function**:

```
Z = Tr(e^{-βH})
```

For a two-level system with Hamiltonian H = diag(0, ε):

```
Z = 1 + e^{-βε}

p₀ = 1/Z          (ground state population)
p₁ = e^{-βε}/Z    (excited state population)
```

**Thermodynamic quantities:**

- **Mean energy**: ⟨H⟩ = ε·p₁ = ε/(1 + e^{βε})
- **Heat capacity**: C = (βε)² · e^{-βε} / (1 + e^{-βε})²
- **Entropy**: S = -p₀ ln p₀ - p₁ ln p₁
- **Free energy**: F = -β⁻¹ ln Z

The heat capacity exhibits the **Schottky anomaly** — a peak at βε ≈ 2.4, characteristic
of two-level systems.

### Quantum Work {#quantum-work-theory}

In quantum mechanics, work is not an observable (a fundamental result by Allahverdyan, 2014).
Instead, we use the **two-point measurement (TPM) scheme**:

1. Measure Hᵢ at time t₀ → obtain eigenvalue Eₘ
2. Evolve the system under some protocol
3. Measure H_f at time t₁ → obtain eigenvalue Eₙ
4. Work: W = Eₙ - Eₘ

The work distribution is:

```
P(W) = Σₘₙ δ(W - (Eₙ - Eₘ)) pₘ P(n|m)
```

where pₘ is the initial thermal population and P(n|m) is the transition probability.

**The Jarzynski equality** (Jarzynski, 1997; extended to quantum by Mukamel, 2003):

```
⟨e^{-βW}⟩ = Z_f / Z_i
```

This remarkable identity holds for *any* driving protocol, far from equilibrium.
It connects a nonequilibrium fluctuation average to an equilibrium free energy difference.

**The Crooks fluctuation theorem**:

```
P_fwd(W) / P_rev(-W) = e^{β(W - ΔF)}
```

This relates the forward and reverse work distributions, with ΔF = F_f - Fᵢ.

### Von Neumann Entropy {#von-neumann-entropy-theory}

The quantum generalization of Shannon entropy:

```
S(ρ) = -Tr(ρ ln ρ) = -Σᵢ λᵢ ln λᵢ
```

where λᵢ are the eigenvalues of ρ. This is the **von Neumann entropy**, introduced by
John von Neumann in his 1932 book *Mathematische Grundlagen der Quantenmechanik*.

**Properties:**
- S(ρ) ≥ 0 (non-negativity)
- S(ρ) = 0 iff ρ is pure (vanishing for pure states)
- S(ρ) ≤ ln d (maximum for maximally mixed state in d dimensions)
- Subadditivity: S(ρ_AB) ≤ S(ρ_A) + S(ρ_B)
- Concave: S(λρ₁ + (1-λ)ρ₂) ≥ λS(ρ₁) + (1-λ)S(ρ₂)

**Conditional entropy**: S(A|B) = S(ρ_AB) - S(ρ_B)

This can be **negative** for entangled states, unlike the classical case!
A negative conditional entropy indicates quantum entanglement.

**Quantum mutual information**: I(A:B) = S(A) + S(B) - S(AB)

Measures the total (classical + quantum) correlations between A and B.
For a Bell state, I(A:B) = 2 ln 2.

**Entanglement entropy**: For a pure bipartite state |ψ_AB⟩:

```
E(|ψ_AB⟩) = S(ρ_A) = S(Tr_B(|ψ_AB⟩⟨ψ_AB|))
```

This equals zero for product states and ln 2 for maximally entangled states.

**Rényi entropy** of order α:

```
S_α(ρ) = (1/(1-α)) ln(Σᵢ λᵢ^α)
```

Recovers von Neumann entropy as α → 1 and the purity as α → 2 (up to a log).

### Resource Theories {#resource-theories-theory}

A **resource theory** is defined by:
1. **Free states**: states that can be created at no cost
2. **Free operations**: transformations that cannot create resources

| Resource Theory | Free States | Free Operations | Resource |
|----------------|-------------|-----------------|----------|
| Entanglement | Separable | LOCC | Entanglement |
| Athermality | Thermal/Gibbs | Thermal operations | Nonequilibrium |
| Coherence | Diagonal | Incoherent operations | Quantum coherence |

**Resource monotones** are functions that never increase under free operations:

- **Entanglement**: entanglement entropy, concurrence, entanglement of formation
- **Athermality**: relative entropy of athermality S(ρ ‖ ρ_th)
- **Coherence**: L1-norm of coherence, relative entropy of coherence

**State conversion** (Nielsen's theorem for entanglement):
ρ → σ by LOCC iff λ_ρ majorizes λ_σ (λ_ρ ≻ λ_σ).

---

## Architecture

```
                    ┌─────────────────────────────────────┐
                    │     quantum-thermo library           │
                    └──────────────┬──────────────────────┘
                                   │
        ┌──────────┬───────────────┼───────────────┬──────────┐
        │          │               │               │          │
   ┌────▼────┐ ┌───▼───┐   ┌──────▼──────┐  ┌─────▼────┐ ┌──▼───┐
   │  qubit  │ │channel│   │   thermal   │  │   work   │ │entropy│
   │         │ │       │   │             │  │          │ │quantum│
   │ Qubit   │ │Kraus  │   │ Thermal     │  │Work      │ │       │
   │ State   │ │ops    │   │ State       │  │Distrib.  │ │Von    │
   │         │ │       │   │             │  │          │ │Neumann│
   │ Density │ │BitFlip│   │ Partition   │  │Jarzynski │ │       │
   │ Matrix  │ │Phase  │   │ Function    │  │Equality  │ │Mutual │
   │         │ │Depol. │   │             │  │          │ │Info   │
   │ Pauli   │ │AmpDamp│   │ Heat        │  │Crooks   │ │       │
   │ Bloch   │ │       │   │ Capacity    │  │Fluct.   │ │Entang.│
   └─────────┘ └───────┘   └─────────────┘  └──────────┘ └───────┘
        │          │               │               │          │
        └──────────┴───────────────┴───────────────┴──────────┘
                                   │
                    ┌──────────────▼──────────────────────┐
                    │          resource module             │
                    │  Entanglement · Athermality ·       │
                    │  Coherence · LOCC · Conversion       │
                    └─────────────────────────────────────┘
```

### Data Flow

```
QubitState ───► DensityMatrix ───► QuantumChannel.apply()
                                       │
                                       ▼
                                  Transformed State
                                       │
                    ┌──────────────────┤
                    │                  │
              ThermalState      WorkDistribution
                    │                  │
                    ▼                  ▼
              Entropy/CV          Jarzynski Equality
                    │
                    ▼
            Resource Monotones
```

---

## Module Reference

| Module | Description | Key Types | Key Functions |
|--------|------------|-----------|---------------|
| `qubit` | Qubit states and density matrices | `QubitState`, `DensityMatrix` | `bloch_coords()`, `purity()`, `eigenvalues()` |
| `channel` | Quantum channels (CPTP maps) | `QuantumChannel` | `bit_flip()`, `depolarizing()`, `amplitude_damping()` |
| `thermal` | Gibbs states and thermodynamics | `ThermalState` | `partition_function()`, `heat_capacity()`, `free_energy()` |
| `work` | Work distributions and fluctuation theorems | `WorkDistribution`, `WorkTrajectory` | `jarzynski_average()`, `mean_work()` |
| `entropy_quantum` | Von Neumann entropy and correlations | `BipartiteState` | `von_neumann_entropy()`, `mutual_information()` |
| `resource` | Quantum resource theories | `ResourceTheoryType`, `ConversionResult` | `resource_monotones()`, `can_convert()` |

---

## Examples

### Example 1: Bloch Sphere Exploration

Creating and measuring qubit states on the Bloch sphere.

```rust
use quantum_thermo::qubit::{QubitState, DensityMatrix, PAULI_X, PAULI_Z};

fn main() {
    // Standard basis states
    let zero = QubitState::zero();   // |0⟩ = north pole
    let one = QubitState::one();     // |1⟩ = south pole
    let plus = QubitState::plus();   // |+⟩ = equator, x-direction
    let minus = QubitState::minus(); // |−⟩ = equator, -x direction

    // Check Bloch sphere coordinates
    let (theta_0, _) = zero.bloch_coords();
    println!("|0⟩ → theta = {:.4} (north pole)", theta_0);

    let (theta_1, _) = one.bloch_coords();
    println!("|1⟩ → theta = {:.4} (south pole)", theta_1);

    // Convert to density matrices
    let rho_plus = plus.to_density_matrix();
    println!("|+⟩ density matrix:");
    println!("  [{:.4}, {:.4}]", rho_plus.data[0][0], rho_plus.data[0][1]);
    println!("  [{:.4}, {:.4}]", rho_plus.data[1][0], rho_plus.data[1][1]);

    // Check purity
    println!("Purity of |+⟩: {:.6} (pure = 1.0)", rho_plus.purity());

    // Create a mixed state: 50/50 mixture of |0⟩ and |1⟩
    let mixed = DensityMatrix::from_mixture(0.5, &zero, &one);
    println!("Mixed state purity: {:.6} (maximally mixed = 0.5)", mixed.purity());

    // Eigenvalues tell us the spectrum
    let (e1, e2) = mixed.eigenvalues();
    println!("Eigenvalues: [{:.4}, {:.4}]", e1, e2);

    // Apply Pauli X (bit flip)
    let flipped = zero.to_density_matrix().apply_matrix(&PAULI_X);
    println!("X|0⟩⟨0|X: [{:.4}, {:.4}]", flipped.data[0][0], flipped.data[1][1]);

    // Create from Bloch angles
    let theta = std::f64::consts::FRAC_PI_4; // 45 degrees
    let state = QubitState::from_bloch(theta, 0.0);
    println!("State at θ=π/4: alpha={:.4}, beta={:.4}", state.alpha, state.beta);
    println!("  P(|0⟩) = {:.4}, P(|1⟩) = {:.4}", state.prob_zero(), state.prob_one());
}
```

### Example 2: Quantum Channels and Fidelity

Applying quantum channels and computing fidelity between states.

```rust
use quantum_thermo::qubit::QubitState;
use quantum_thermo::channel::{QuantumChannel, fidelity, trace_distance};

fn main() {
    // Start with |+⟩ state
    let psi = QubitState::plus();
    let rho = psi.to_density_matrix();
    println!("Initial state |+⟩:");
    println!("  Purity: {:.6}", rho.purity());

    // Apply bit-flip channel with increasing noise
    println!("\n--- Bit Flip Channel ---");
    for p in [0.0, 0.1, 0.25, 0.5, 1.0] {
        let ch = QuantumChannel::bit_flip(p);
        let rho_out = ch.apply(&rho);
        let f = fidelity(&rho, &rho_out);
        println!("p={:.2}: fidelity={:.4}, purity={:.4}",
                 p, f, rho_out.purity());
    }

    // Amplitude damping: models energy dissipation
    println!("\n--- Amplitude Damping ---");
    let excited = QubitState::one().to_density_matrix();
    for gamma in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let ch = QuantumChannel::amplitude_damping(gamma);
        let rho_out = ch.apply(&excited);
        println!("γ={:.2}: ⟨1|ρ|1⟩={:.4} (excited population)",
                 gamma, rho_out.data[1][1]);
    }

    // Depolarizing channel
    println!("\n--- Depolarizing Channel ---");
    let ch = QuantumChannel::depolarizing(0.5);
    let rho_depol = ch.apply(&rho);
    println!("p=0.5: fidelity={:.4}", fidelity(&rho, &rho_depol));
    println!("  Approaches I/2 for large p");

    // Trace distance between |0⟩ and |1⟩
    let rho0 = QubitState::zero().to_density_matrix();
    let rho1 = QubitState::one().to_density_matrix();
    let td = trace_distance(&rho0, &rho1);
    println!("\nTrace distance |0⟩↔|1⟩: {:.4} (maximum = 1.0)", td);

    // Verify channel completeness
    let ch = QuantumChannel::depolarizing(0.3);
    println!("Depolarizing(0.3) completeness: {}", ch.verify_completeness());
}
```

### Example 3: Thermal States at Various Temperatures

Computing thermal states and thermodynamic quantities.

```rust
use quantum_thermo::thermal::ThermalState;

fn main() {
    let energy_gap = 1.0; // Energy gap ε = 1.0 (natural units)

    println!("Two-level system with energy gap ε = 1.0");
    println!("{:-<60}");
    println!("{:>8} {:>8} {:>8} {:>8} {:>10}",
             "T", "p0", "p1", "<E>", "S");
    println!("{:-<60}");

    for beta in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0] {
        let t = 1.0 / beta;
        let ts = ThermalState::new(energy_gap, beta);
        println!("{:>8.3} {:>8.4} {:>8.4} {:>8.4} {:>10.4}",
                 t,
                 ts.ground_population(),
                 ts.excited_population(),
                 ts.mean_energy(),
                 ts.entropy());
    }

    // Detailed analysis at β = 1
    println!("\n--- Detailed Analysis at β = 1 ---");
    let ts = ThermalState::new(1.0, 1.0);
    println!("Partition function Z = {:.6}", ts.partition_function());
    println!("Free energy F = {:.6}", ts.free_energy());
    println!("Internal energy U = {:.6}", ts.internal_energy());
    println!("Heat capacity C = {:.6}", ts.heat_capacity());
    println!("Entropy S = {:.6}", ts.entropy());

    // Heat capacity scan: Schottky anomaly
    println!("\n--- Heat Capacity (Schottky Anomaly) ---");
    let mut max_c = 0.0;
    let mut max_t = 0.0;
    for i in 1..100 {
        let t = 0.01 * i as f64;
        let ts = ThermalState::from_temperature(1.0, t);
        let c = ts.heat_capacity();
        if c > max_c {
            max_c = c;
            max_t = t;
        }
    }
    println!("Peak heat capacity: C = {:.6} at T = {:.4}", max_c, max_t);
    println!("(Theoretical peak at βε ≈ 2.4, T ≈ 0.42)");
}
```

### Example 4: Von Neumann Entropy and Entanglement

Computing entropic quantities and quantifying entanglement.

```rust
use quantum_thermo::qubit::{QubitState, DensityMatrix};
use quantum_thermo::entropy_quantum::{
    von_neumann_entropy, mutual_information, entanglement_entropy,
    renyi_entropy, BipartiteState,
};
use quantum_thermo::resource::{entanglement_monotones, is_separable, ResourceTheoryType};

fn main() {
    // Von Neumann entropy for different states
    println!("--- Von Neumann Entropy ---");
    let rho_pure = QubitState::plus().to_density_matrix();
    println!("Pure |+⟩: S = {:.6} (should be 0)",
             von_neumann_entropy(&rho_pure));

    let rho_mixed = DensityMatrix::maximally_mixed();
    println!("Max mixed: S = {:.6} (should be ln 2 ≈ 0.693)",
             von_neumann_entropy(&rho_mixed));

    // Rényi entropy family
    println!("\n--- Rényi Entropy (maximally mixed) ---");
    let eigs = [0.5, 0.5];
    for alpha in [0.5, 1.0, 2.0, 5.0, 100.0] {
        let s = renyi_entropy(&eigs, alpha);
        println!("α={:>6.1}: S_α = {:.6}", alpha, s);
    }

    // Bipartite entanglement
    println!("\n--- Entanglement Analysis ---");
    let bell = BipartiteState::bell_state();
    let product = BipartiteState::product_00();

    let ee_bell = entanglement_entropy(&bell);
    let ee_prod = entanglement_entropy(&product);
    println!("Bell |Φ+⟩ entanglement: {:.6} (should be ln 2)", ee_bell);
    println!("Product |00⟩ entanglement: {:.6} (should be 0)", ee_prod);

    let mi_bell = mutual_information(&bell);
    let mi_prod = mutual_information(&product);
    println!("Bell mutual information: {:.6}", mi_bell);
    println!("Product mutual information: {:.6}", mi_prod);

    // Separability test
    println!("\n--- Separability ---");
    println!("Bell state separable? {}", is_separable(&bell));
    println!("Product state separable? {}", is_separable(&product));

    // Entanglement monotones for Bell state
    let monos = entanglement_monotones(&bell);
    println!("\n--- Bell State Monotones ---");
    for m in &monos {
        println!("{}: {:.6}", m.name, m.value);
    }

    // Resource monotones for single-qubit coherence
    let plus = QubitState::plus().to_density_matrix();
    let monos = quantum_thermo::resource::resource_monotones(
        &plus, ResourceTheoryType::Coherence);
    println!("\n--- Coherence of |+⟩ ---");
    for m in &monos {
        println!("{}: {:.6}", m.name, m.value);
    }
}
```

---

## Performance Analysis

### Time Complexity

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| Density matrix creation | O(1) | Fixed 2×2 matrix |
| Purity computation | O(1) | 4 multiplications + 2 additions |
| Eigenvalues (2×2) | O(1) | Closed-form: discriminant formula |
| Channel application (k Kraus ops) | O(k) | k matrix multiplications |
| Von Neumann entropy | O(1) | 2 eigenvalues + 2 log evaluations |
| Power method (4×4, n iterations) | O(n · 16) | n iterations × 16 multiply-adds |
| Partial trace | O(1) | 8 element accesses |
| Bloch coordinates | O(1) | 1 acos call |

### Space Complexity

All types use fixed-size allocations:
- `DensityMatrix`: 4 × 8 bytes = 32 bytes
- `QubitState`: 2 × 8 bytes = 16 bytes
- `BipartiteState`: 16 × 8 bytes = 128 bytes
- `QuantumChannel`: up to 4 × 32 bytes = 128 bytes (Kraus ops)

### Numerical Stability

- **Eigenvalue computation**: For 2×2 matrices, we use the closed-form
  discriminant formula which is numerically stable for well-conditioned matrices.
- **Entropy computation**: Eigenvalues below 1e-15 are treated as zero to avoid
  log(0) = -∞. This is physically correct since zero eigenvalues contribute
  nothing to the entropy (lim_{x→0} x ln x = 0).
- **Kraus completeness**: Verified to 1e-10 tolerance. Small deviations are
  expected from floating-point arithmetic.
- **Power method**: Converges geometrically with rate |λ₂/λ₁|. For degenerate
  eigenvalues, convergence may be slow — the iteration count parameter allows
  the user to trade speed for accuracy.

---

## Design Decisions

### Why Real Matrix Representation?

We use `[[f64; 2]; 2]` (real-valued 2×2 matrices) instead of a complex number library
for several reasons:

1. **Zero dependencies**: No external complex number crate needed. Only `serde` for
   serialization.
2. **Sufficiency for 2-level systems**: Many quantum thermodynamic calculations
   involve diagonal density matrices (thermal states) or real-valued states (|0⟩, |1⟩,
   |+⟩, |−⟩). The full complex structure is unnecessary for these cases.
3. **Performance**: Real arithmetic is ~2× faster than complex arithmetic with the
   same operation count, and has better cache behavior.
4. **Clarity**: Working with real numbers makes the physics more transparent. You can
   see exactly what's happening without complex number bookkeeping.

**Limitation**: States with non-trivial phases (e.g., |+i⟩ = (|0⟩ + i|1⟩)/√2) cannot
be fully represented. The Pauli Y matrix is stored as a zero matrix since its off-diagonal
elements are purely imaginary.

### Why 2-Level Systems Only?

Quantum thermodynamics is most tractable for two-level systems, which are also
the most physically relevant for near-term quantum technologies:

1. **Analytical results**: Most quantities (entropy, fidelity, eigenvalues) have
   closed-form expressions for 2×2 matrices.
2. **Qubit focus**: Quantum computing operates on qubits. A 2-level library
   directly models qubit thermodynamics.
3. **Pedagogical value**: All the key concepts — entanglement, decoherence,
   fluctuation theorems — can be demonstrated with two-level systems.

### Why Iterative Eigenvalue Computation?

For 4×4 matrices (bipartite systems), we use the power method rather than
recursive QR iteration:

1. **Simplicity**: Power method is straightforward to implement and understand.
2. **No external LAPACK dependency**: Keeps the crate self-contained.
3. **Sufficient accuracy**: For the symmetric matrices encountered in quantum
   mechanics, 50-100 iterations give excellent accuracy.
4. **Predictable performance**: Fixed iteration count = predictable runtime.

### Why No `num-complex` or `nalgebra`?

This crate explicitly avoids heavy linear algebra dependencies. For 2×2 and 4×4
matrices, hand-written operations are:

- Faster to compile (no generic linear algebra)
- Easier to audit (all code is visible)
- More portable (no system libraries needed)
- Sufficient for the target use case

---

## Classical vs Quantum Comparison

### Entropy

| Property | Classical (Shannon) | Quantum (von Neumann) |
|----------|-------------------|----------------------|
| Definition | H(p) = -Σ pᵢ log pᵢ | S(ρ) = -Tr(ρ ln ρ) |
| Input | Probability distribution | Density matrix |
| Range | [0, log d] | [0, log d] |
| Zero when | Delta distribution | Pure state |
| Maximum | Uniform distribution | Maximally mixed state |
| Conditioning | H(X|Y) ≥ 0 always | S(A|B) can be **negative** |
| Subadditivity | H(XY) ≤ H(X) + H(Y) | S(ρ_AB) ≤ S(ρ_A) + S(ρ_B) |

### Channels

| Property | Classical Channel | Quantum Channel |
|----------|------------------|-----------------|
| Representation | Stochastic matrix | Kraus operators |
| Input | Probability vector | Density matrix |
| Output | Probability vector | Density matrix |
| Capacity | Mutual information maximized | Holevo quantity |
| Noise model | Bit flip (probabilistic) | Bit flip (coherent + probabilistic) |
| Key difference | Only diagonal elements matter | Off-diagonal (coherence) matters |

### Key Quantum Advantages

1. **Negative conditional entropy**: Enables quantum superdense coding
2. **Quantum mutual information**: Can exceed classical bounds
3. **Coherence as resource**: Off-diagonal elements carry thermodynamic value
4. **Entanglement**: Enables correlations stronger than any classical mechanism

---

## Glossary

**Amplitude damping**: Quantum channel modeling spontaneous emission. An excited
state |1⟩ decays to |0⟩ at rate γ, transferring energy to the environment.

**Bloch sphere**: Geometric representation of a qubit state as a point on or inside
a unit sphere. Pure states are on the surface; mixed states are inside.

**CPTP map**: Completely Positive Trace-Preserving map. The most general physical
transformation of a quantum state. Equivalent to a quantum channel.

**Depolarizing channel**: Noise model that replaces the state with the maximally
mixed state I/2 with probability p, leaving it unchanged with probability 1-p.

**Density matrix**: Hermitian positive semidefinite matrix with unit trace,
describing the statistical state of a quantum system. Generalizes the pure
state |ψ⟩ to mixed states.

**Entanglement entropy**: Von Neumann entropy of the reduced density matrix of
a bipartite pure state. Measures quantum correlations between subsystems.

**Fluctuation theorem**: Symmetry relation between the probability of forward
and reverse trajectories. Implies the second law of thermodynamics as a
statistical average.

**Free energy**: Thermodynamic potential F = ⟨H⟩ - TS. Minimized at thermal
equilibrium.

**Gibbs state**: Thermal equilibrium state ρ = e^{-βH}/Z. Maximizes entropy
at fixed mean energy.

**Jarzynski equality**: ⟨e^{-βW}⟩ = Z_f/Z_i. An exact identity relating
nonequilibrium work fluctuations to equilibrium free energy differences.

**Kraus operator**: Component of the Kraus representation of a quantum channel.
Satisfies Σ_k E_k† E_k = I.

**LOCC**: Local Operations and Classical Communication. The free operations in
entanglement theory. Two parties can act locally and share classical information
but cannot create entanglement.

**Partition function**: Z = Tr(e^{-βH}). Central object in statistical mechanics.
Encodes all thermodynamic information about a system at temperature T.

**Pauli matrices**: The three 2×2 traceless Hermitian matrices X, Y, Z. Together
with the identity, form a basis for all 2×2 Hermitian matrices.

**Purity**: Tr(ρ²). Equals 1 for pure states, less than 1 for mixed states.
Related to the length of the Bloch vector: Tr(ρ²) = (1 + |r|²)/2.

**Resource monotone**: Function that never increases under free operations.
Used to quantify the "value" of a resource state.

**Schottky anomaly**: Peak in the heat capacity of a two-level system at
temperature T ≈ ε/(2.4 k_B). Arises from the thermal activation of the excited state.

**Thermal state**: See Gibbs state.

**Trace distance**: ½|ρ - σ|₁. Metric on density matrices with operational
meaning: maximum probability of distinguishing ρ from σ in a single measurement.

**Two-point measurement**: Work definition scheme using projective measurements
before and after a thermodynamic process. Avoids the problem that quantum work
is not an observable.

**Von Neumann entropy**: S(ρ) = -Tr(ρ ln ρ). Quantum generalization of Shannon
entropy. Measures the uncertainty in the state of a quantum system.

---

## References

1. **von Neumann, J.** (1932). *Mathematische Grundlagen der Quantenmechanik*.
   Springer, Berlin. — The original formulation of quantum statistical mechanics,
   including the density matrix formalism and von Neumann entropy.

2. **Landauer, R.** (1961). "Irreversibility and heat generation in the computing
   process." *IBM Journal of Research and Development*, 5(3), 183–191. — Established
   the fundamental link between information erasure and thermodynamic heat dissipation:
   erasing one bit requires at least k_B T ln 2 of energy.

3. **Bennett, C. H.** (1982). "The thermodynamics of computation — a review."
   *International Journal of Theoretical Physics*, 21(12), 905–940. — Showed that
   computation can be performed reversibly, overturning the belief that computation
   inherently dissipates energy.

4. **Jarzynski, C.** (1997). "Nonequilibrium equality for free energy differences."
   *Physical Review Letters*, 78(14), 2690. — The landmark paper establishing the
   Jarzynski equality, connecting nonequilibrium fluctuations to equilibrium
   thermodynamics.

5. **Scully, M. O.** (2002). "Quantum afterburner: Improving the efficiency of an
   ideal heat engine." *Physical Review Letters*, 88(5), 050602. — Demonstrated that
   quantum coherence can be used to extract more work from a heat engine than
   classical thermodynamics allows.

6. **Goold, J., Huber, M., Riera, A., del Rio, L., & Skrzypczyk, P.** (2016).
   "The role of quantum information in thermodynamics — a topical review."
   *Journal of Physics A: Mathematical and Theoretical*, 49(14), 143001. —
   Comprehensive review of quantum thermodynamics, covering resource theories,
   fluctuation theorems, and quantum thermal machines.

7. **Deffner, S. & Campbell, S.** (2019). *Quantum Thermodynamics: An Introduction
   to the Thermodynamics of Quantum Information*. Morgan & Claypool. — Modern
   textbook covering quantum speed limits, fluctuation theorems, and the thermodynamic
   cost of quantum measurements.

8. **Nielsen, M. A. & Chuang, I. L.** (2010). *Quantum Computation and Quantum
   Information* (10th anniversary edition). Cambridge University Press. — The
   standard reference for quantum information theory, including density matrices,
   quantum channels, and entanglement theory.

9. **Brandão, F. G. S. L., Horodecki, M., Oppenheim, J., Renes, J. M., & Spekkens,
   R. W.** (2013). "Resource theory of quantum states out of thermal equilibrium."
   *Physical Review Letters*, 111(25), 250404. — Formalized the resource theory of
   athermality and its connections to thermodynamic transformations.

10. **Lostaglio, M., Jennings, D., & Rudolph, T.** (2015). "Description of quantum
    coherence in thermodynamic processes requires constraints beyond free energy."
    *Nature Communications*, 6, 6383. — Showed that quantum coherence imposes
    additional thermodynamic constraints beyond what free energy alone predicts.

---

## License

MIT
