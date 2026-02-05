import pennylane as qml
from pennylane import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

print("PennyLane version:", qml.__version__)
print("Environment setup successful!")

symbols = ["H", "H"]
coordinates = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.4])

H_h2, qubits = qml.qchem.molecular_hamiltonian(symbols,coordinates,charge=0,mult=1,basis='sto-3g')

print(f"Molecular Hamiltonian for H₂:")
print(f"Number of qubits required: {qubits}")

if hasattr(H_h2,'ops'):
    num_terms = len(H_h2.ops)
    coeffs = H_h2.coeffs[:5]
    ops = H_h2.ops[:5]
else:
    if hasattr(H_h2, 'operands'):
        num_terms = len(H_h2.operands)
        # Get first 5 terms
        terms_to_show = H_h2.operands[:5]
    else:
        # Alternative: get pauli representation
        num_terms = len(H_h2.pauli_rep) if hasattr(H_h2, 'pauli_rep') else "N/A"
        terms_to_show = list(H_h2.pauli_rep.items())[:5] if hasattr(H_h2, 'pauli_rep') else []

print(f"Number of terms in Hamiltonian: {num_terms}")
print(f"\nHamiltonian (first 5 terms):")

if hasattr(H_h2, 'ops'):
    for i, (coeff, op) in enumerate(zip(coeffs, ops)):
        print(f"  {coeff.real:.6f} * {op}")
elif hasattr(H_h2, 'operands'):
    for i, term in enumerate(terms_to_show):
        print(f"  {term}")
else:
    for pw, coeff in terms_to_show:
        print(f"  {coeff:.6f} * {pw}")

dev = qml.device('default.qubit', wires=qubits)

def hardware_efficient_ansatz(params, wires):
    n_qubits = wires
    n_layers = len(params)

    for layer in range(n_layers):
        for i in range(n_qubits):
            qml.RY(params[layer][i], wires=i)
            qml.RZ(params[layer][i + n_qubits], wires=i)

        # Full entanglement: circular CNOT chain
        for i in range(n_qubits - 1):
            qml.CNOT(wires=[i, i + 1])
        if n_qubits > 2:
            qml.CNOT(wires=[n_qubits - 1, 0])

@qml.qnode(dev)
def vqe_circuit(params, hamiltonian):
    # Hartree-Fock state for H2: |1100>
    qml.PauliX(wires=0)
    qml.PauliX(wires=1)
    hardware_efficient_ansatz(params, qubits)
    return qml.expval(hamiltonian)

n_layers = 4
params_shape = (n_layers, 2 * qubits)
np.random.seed(42)
init_params = np.random.random(params_shape) * 0.1

print("VQE Circuit Structure:")
print(qml.draw(vqe_circuit)(init_params, H_h2))

def cost_function(params_flat):
    params = params_flat.reshape(params_shape)
    return vqe_circuit(params, H_h2)

init_params_flat = init_params.flatten()
energy_history = []

def callback(params):
    energy = cost_function(params)
    energy_history.append(energy)
    if len(energy_history) % 10 == 0:
        print(f"Iteration {len(energy_history)}: Energy = {energy:.8f} Ha")

print("Starting VQE optimization...\n")
print(f"Initial energy: {cost_function(init_params_flat):.8f} Hartree\n")

result = minimize(
    cost_function,
    init_params_flat,
    method='L-BFGS-B',
    jac='3-point',
    callback=callback,
    options={'maxiter': 300, 'ftol': 1e-12}
)

vqe_energy = result.fun
print(f"\n{'='*50}")
print(f"VQE Optimization Complete!")
print(f"Final Ground State Energy: {vqe_energy:.8f} Hartree")
print(f"{'='*50}")

import numpy as np_regular
mol = qml.qchem.Molecule(symbols, np_regular.array(coordinates.reshape(-1, 3)), charge=0, mult=1, basis_name='sto-3g')
hf_energy = qml.qchem.scf(mol)()[0][0]

exact_energy = -1.137283

vqe_error = abs(vqe_energy - exact_energy)
hf_error = abs(hf_energy - exact_energy)

print("\n" + "="*60)
print("COMPARISON WITH CLASSICAL METHODS")
print("="*60)
print(f"\nExact Ground State Energy (FCI): {exact_energy:.8f} Ha")
print(f"Hartree-Fock Energy (Classical):  {hf_energy:.8f} Ha")
print(f"VQE Energy (Quantum):             {vqe_energy:.8f} Ha")
print(f"\nAbsolute Errors:")
print(f"  Hartree-Fock Error: {hf_error:.8f} Ha ({hf_error*627.5:.4f} kcal/mol)")
print(f"  VQE Error:          {vqe_error:.8f} Ha ({vqe_error*627.5:.4f} kcal/mol)")
print(f"\nImprovement over HF: {((hf_error - vqe_error)/hf_error * 100):.2f}%")

if vqe_error < hf_error:
    print("\n✓ VQE achieved better accuracy than Hartree-Fock!")
else:
    print("\n⚠ VQE needs more optimization layers or iterations")

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(energy_history, 'b-', linewidth=2, label='VQE Energy')
plt.axhline(y=exact_energy, color='r', linestyle='--', linewidth=2, label='Exact Energy')
plt.axhline(y=hf_energy, color='g', linestyle='--', linewidth=2, label='Hartree-Fock')
plt.xlabel('Optimization Iteration', fontsize=12)
plt.ylabel('Energy (Hartree)', fontsize=12)
plt.title('VQE Energy Convergence', fontsize=14, fontweight='bold')
plt.legend(fontsize=10)
plt.grid(True, alpha=0.3)

plt.subplot(1, 2, 2)
errors = [abs(e - exact_energy) for e in energy_history]
plt.semilogy(errors, 'b-', linewidth=2)
plt.xlabel('Optimization Iteration', fontsize=12)
plt.ylabel('Absolute Error (Hartree)', fontsize=12)
plt.title('VQE Error Convergence (Log Scale)', fontsize=14, fontweight='bold')
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('vqe_convergence.png', dpi=300, bbox_inches='tight')
plt.show()

print("\nConvergence plot saved!")

bond_lengths = np.linspace(0.5, 4.5, 15)
vqe_energies = []
hf_energies = []

print("Scanning potential energy surface...\n")

# Use warm-starting: previous bond length's optimal params as initial guess for next
current_params = result.x.copy()

for bond_length in bond_lengths:
    coords = np.array([0.0, 0.0, 0.0, 0.0, 0.0, bond_length])

    H, _ = qml.qchem.molecular_hamiltonian(symbols, coords, basis='sto-3g')

    # Re-optimize VQE for each bond length (warm-started from previous)
    def cost_fn_scan(params_flat):
        params = params_flat.reshape(params_shape)
        return vqe_circuit(params, H)

    res_scan = minimize(
        cost_fn_scan,
        current_params,
        method='L-BFGS-B',
        jac='3-point',
        options={'maxiter': 200, 'ftol': 1e-12}
    )
    energy_vqe = res_scan.fun
    current_params = res_scan.x.copy()  # warm-start next iteration
    vqe_energies.append(energy_vqe)

    mol_scan = qml.qchem.Molecule(symbols, np_regular.array(coords.reshape(-1, 3)), basis_name='sto-3g')
    energy_hf = qml.qchem.scf(mol_scan)()[0][0]
    hf_energies.append(energy_hf)

    print(f"Bond length: {bond_length:.2f} Bohr ({bond_length*0.529:.2f} Å) | VQE: {energy_vqe:.6f} Ha | HF: {energy_hf:.6f} Ha")

plt.figure(figsize=(10, 6))
plt.plot(bond_lengths * 0.529, vqe_energies, 'bo-', linewidth=2, markersize=8, label='VQE')
plt.plot(bond_lengths * 0.529, hf_energies, 'gs--', linewidth=2, markersize=8, label='Hartree-Fock')
plt.xlabel('H-H Bond Length (Angstroms)', fontsize=12)
plt.ylabel('Energy (Hartree)', fontsize=12)
plt.title('H₂ Potential Energy Surface', fontsize=14, fontweight='bold')
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.savefig('h2_potential_energy_surface.png', dpi=300, bbox_inches='tight')
plt.show()

print(f"\nOptimal bond length: {bond_lengths[np.argmin(vqe_energies)] * 0.529:.3f} Å")
print(f"Minimum energy: {min(vqe_energies):.8f} Ha")