# VQE for Molecular Ground State Energy Calculation

**Name:** SK SAIF IBNA EZHAR ARKO  
**Objective:** Implement Variational Quantum Eigensolver (VQE) to calculate ground state energy of small molecules

## Overview

This project implements the Variational Quantum Eigensolver (VQE) algorithm using PennyLane to calculate the ground state energy of the hydrogen molecule (H₂). VQE is a hybrid quantum-classical algorithm particularly suitable for near-term quantum devices (NISQ era) that can achieve chemical accuracy in molecular energy calculations.

## Key Concepts

- **Molecular Hamiltonian**: Quantum operator describing molecular energy
- **Ansatz**: Parameterized quantum circuit that prepares trial wavefunctions
- **Variational Principle**: Ground state has the lowest energy expectation value
- **Hybrid Optimization**: Quantum circuit evaluation + classical parameter optimization

## Environment Setup

The code requires PennyLane and its dependencies. Install with:
```bash
pip install pennylane matplotlib scipy
```

## Implementation Structure

### 1. Molecular Hamiltonian Construction
- Defines H₂ molecule with bond length of 1.4 Bohr
- Uses STO-3G basis set
- Generates molecular Hamiltonian using `qml.qchem.molecular_hamiltonian()`
- Displays first 5 terms of the Hamiltonian for verification

### 2. VQE Circuit Design

#### Hardware Efficient Ansatz (HEA)
The ansatz structure consists of 4 layers, each containing:
- **RY rotation gates** on all qubits (parameterized)
- **RZ rotation gates** on all qubits (parameterized)
- **Circular CNOT chain** for entanglement (including final-to-first qubit connection)

#### VQE Quantum Circuit
The complete circuit includes:
- **Initial state preparation**: Hartree-Fock state for H₂ (|1100⟩ - both electrons in bonding orbital)
- **Parameterized ansatz**: Hardware-efficient ansatz with trainable parameters
- **Energy measurement**: Expectation value of molecular Hamiltonian

### 3. Classical Optimization
- Uses **L-BFGS-B** optimizer from SciPy
- Employs 3-point numerical gradient approximation
- Maximum 300 iterations with tolerance of 1e-12
- Callback function tracks energy convergence history

### 4. Classical Method Comparison
Compares VQE results with:
- **Hartree-Fock energy**: Mean-field approximation (classical reference)
- **Exact ground state energy**: Full Configuration Interaction (FCI) at -1.137283 Ha
- Calculates absolute errors and improvement percentage

### 5. Visualization
Generates two plots showing:
- **Energy convergence**: VQE energy vs. optimization iterations with HF and exact reference lines
- **Error convergence**: Absolute error on logarithmic scale
- Saves as `vqe_convergence.png` at 300 DPI

### 6. Potential Energy Surface Scan
- Scans 15 bond lengths from 0.5 to 4.5 Bohr
- Uses **warm-starting**: optimal parameters from previous bond length as initial guess
- Compares VQE vs. Hartree-Fock across the potential energy surface
- Maximum 200 iterations per bond length
- Generates `h2_potential_energy_surface.png`
- Identifies optimal bond length and minimum energy

## Results and Key Findings

### What the Implementation Demonstrates:

1. **VQE Algorithm End-to-End**
   - Constructs molecular Hamiltonian for H₂ with 4 qubits
   - Implements hardware-efficient ansatz with 4 layers
   - Optimizes circuit parameters using L-BFGS-B
   - Achieves ground state energy calculation with chemical accuracy

2. **Quantum Advantage Over Classical Methods**
   - VQE captures electron correlation effects beyond Hartree-Fock
   - Shows percentage improvement over mean-field approximation
   - Demonstrates convergence to exact FCI energy
   - Provides validation against classical quantum chemistry

3. **Molecular Binding Behavior**
   - Potential energy surface reveals H-H bond dissociation
   - Identifies equilibrium bond length (~0.74 Å for H₂)
   - Compares quantum vs. classical predictions across bond lengths
   - Shows characteristic binding energy curve

4. **Practical Quantum Computing**
   - Uses `default.qubit` simulator (4 qubits required)
   - Prints circuit structure for transparency
   - Tracks optimization progress with callbacks
   - Generates publication-quality visualizations

3. **Drug Discovery Relevance**
   - Accurate energy calculations crucial for molecular stability
   - Understanding binding energies in drug-protein interactions
   - Predicting chemical reaction pathways

## Limitations and Considerations

- **Qubit Requirements**: Current implementation uses 4 qubits for H₂; larger molecules scale exponentially
- **Basis Set**: STO-3G is minimal; larger basis sets improve accuracy but increase qubits needed
- **Optimizer**: L-BFGS-B can get stuck in local minima; multiple random initializations may be needed
- **Noise**: Simulation is noise-free; real quantum hardware (NISQ devices) will have lower accuracy
- **Classical Overhead**: Numerical gradients require multiple circuit evaluations per iteration
- **Ansatz Design**: Hardware-efficient ansatz is general but may not be optimal for all molecules

## Next Steps for Extension

- **Larger Molecules**: Extend to H₂O, NH₃, or small organic molecules
- **Quantum Hardware**: Test on IBM Quantum, IonQ, or other cloud quantum computers
- **Ansatz Optimization**: Try UCCSD (Unitary Coupled Cluster) or adaptive VQE
- **Integration with Drug Discovery**: Combine with QAOA for molecular docking simulations
- **Excited States**: Implement VQE variants for excited state calculations
- **Error Mitigation**: Add noise models and error mitigation techniques for hardware deployment

## Output Files

The implementation generates two visualization files:

1. **`vqe_convergence.png`**
   - Left plot: Energy convergence showing VQE trajectory with HF and exact energy reference lines
   - Right plot: Absolute error on logarithmic scale tracking convergence to exact solution

2. **`h2_potential_energy_surface.png`**
   - Comparison of VQE vs. Hartree-Fock energies across 15 bond lengths (0.5-4.5 Bohr)
   - Demonstrates dissociation behavior and identifies equilibrium bond length
   - Shows quantum advantage in capturing correlation at stretched geometries

## Running the Code

```bash
python vqe_molecular_energy.py
```

Expected output includes:
- PennyLane version confirmation
- Hamiltonian details (number of qubits and first 5 terms)
- Circuit structure visualization
- Optimization progress (energy every 10 iterations)
- Final comparison table with HF, VQE, and exact energies
- Error analysis and improvement percentage
- Bond length scan progress
- Optimal bond length and minimum energy values
