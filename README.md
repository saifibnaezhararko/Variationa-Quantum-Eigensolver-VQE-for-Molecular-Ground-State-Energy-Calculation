# VQE for Molecular Ground State Energy Calculation

**Name:** SK SAIF IBNA EZHAR ARKO  
**Objective:** Implement Variational Quantum Eigensolver (VQE) to calculate ground state energy of small molecules

## Overview

The Variational Quantum Eigensolver (VQE) is a hybrid quantum-classical algorithm used to find the ground state energy of molecular systems. It's particularly useful for drug discovery as accurate energy calculations are crucial for understanding molecular stability and interactions.

## Key Concepts

- **Molecular Hamiltonian**: Quantum operator describing molecular energy
- **Ansatz**: Parameterized quantum circuit that prepares trial wavefunctions
- **Variational Principle**: Ground state has the lowest energy expectation value
- **Hybrid Optimization**: Quantum circuit evaluation + classical parameter optimization

## Environment Setup

The code uses PennyLane for quantum computing simulations. Install with:
```bash
pip install pennylane
```

## Implementation Structure

### 1. Environment Setup - Basic Quantum Chemistry
Creates a simple molecular Hamiltonian for H₂ (hydrogen molecule) to verify the environment.

### 2. VQE Implementation

#### Hardware Efficient Ansatz (HEA)
The ansatz structure includes:
- Layer of RY rotations on all qubits
- Layer of RZ rotations on all qubits
- Entangling CNOT gates in a linear chain
- Multiple layers for depth

#### VQE Circuit
Returns energy expectation value with:
- Initial state preparation (Hartree-Fock for H₂: |01⟩ - one electron in each orbital)
- Parameterized ansatz application
- Energy expectation measurement

### 3. Parameter Optimization
Uses classical optimization (COBYLA method) to minimize the energy expectation value.

### 4. Classical Method Comparison
Compares VQE results with:
- Hartree-Fock energy (classical reference)
- Exact ground state energy from Full Configuration Interaction (FCI)

### 5. Visualization
- Energy convergence during optimization
- Error convergence on logarithmic scale

### 6. Potential Energy Surface Scan
Calculates energy at different H-H bond lengths to demonstrate molecular binding behavior.

## Results and Key Findings

### What We Demonstrated:

1. **VQE Algorithm Implementation**
   - Built molecular Hamiltonian for H₂
   - Designed hardware-efficient ansatz
   - Optimized quantum circuit parameters
   - Achieved ground state energy calculation

2. **Quantum Advantage**
   - VQE can achieve accuracy beyond Hartree-Fock
   - Captures electron correlation effects
   - Scalable to larger molecules (with more qubits)

3. **Drug Discovery Relevance**
   - Accurate energy calculations crucial for molecular stability
   - Understanding binding energies in drug-protein interactions
   - Predicting chemical reaction pathways

## Limitations

- Current quantum hardware has limited qubits (restricts molecule size)
- Noise in NISQ devices affects accuracy
- Classical optimization can get stuck in local minima
- Basis set limitations in current implementation

## Next Steps

- Apply to larger drug-like molecules
- Test on actual quantum hardware
- Integrate with molecular docking (QAOA - next task)

## Output Files

The code generates:
- `vqe_convergence.png` - Energy and error convergence plots
- `h2_potential_energy_surface.png` - H₂ binding energy curve
