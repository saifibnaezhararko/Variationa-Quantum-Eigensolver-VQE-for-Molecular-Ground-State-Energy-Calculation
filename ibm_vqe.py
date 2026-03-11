# IBM Quantum VQE for H2 Molecular Ground State Energy
#
# USAGE:
#   1. Set IBM_TOKEN below (or export IBM_TOKEN env variable) to your IBM Quantum API token.
#      Tokens are available at: https://quantum.ibm.com/
#   2. (Optional) Set IBM_INSTANCE to your hub/group/project string.
#      Leave as None to auto-discover your first available instance.
#   3. Run:  python ibm_vqe.py
#
# CHANNEL NOTES:
#   - "ibm_quantum" : IBM Quantum Network access (token + hub/group/project instance)
#   - "ibm_cloud"   : IBM Cloud access (API key + CRN service instance)
#
# If no IBM connection can be established the script falls back to a local
# PennyLane simulation so you can still see VQE results.

import os
import warnings
from typing import Optional
import numpy as np_std                  # standard NumPy (not wrapped by PennyLane)
import pennylane as qml
from pennylane import numpy as np       # PennyLane's differentiable numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize

warnings.filterwarnings('ignore')

print("PennyLane version:", qml.__version__)
print("IBM Quantum VQE for H2 Ground State Energy")

# ============================================================
# CONFIG – edit these values or set the matching env variables
# ============================================================
IBM_TOKEN = os.environ.get("IBM_TOKEN", "YOUR_IBM_TOKEN_HERE")

# Valid channel names for qiskit-ibm-runtime:
#   "ibm_quantum"  – IBM Quantum Network (hub/group/project)
#   "ibm_cloud"    – IBM Cloud (CRN-based)
# Do NOT use "ibm_quantum_platform" – that name is invalid.
IBM_CHANNEL = os.environ.get("IBM_CHANNEL", "ibm_quantum")

# hub/group/project instance string for the "ibm_quantum" channel.
# Set to None to auto-discover your first available instance.
# Example: "ibm-q/open/main" (if your account still has access to this instance)
IBM_INSTANCE = os.environ.get("IBM_INSTANCE", None)

# Backend name.  Use "ibmq_qasm_simulator" for cloud simulation or a real
# device name such as "ibm_brisbane".  None = least-busy device.
IBM_BACKEND = os.environ.get("IBM_BACKEND", "ibmq_qasm_simulator")

# When True the script continues with a local simulator if IBM is unavailable.
FALLBACK_TO_LOCAL = True
# ============================================================

# ------------------------------------------------------------------
# Token display helper
# ------------------------------------------------------------------
def _mask_token(token: str) -> str:
    """Return a partially masked token safe to print."""
    if not token or len(token) < 16 or token == "YOUR_IBM_TOKEN_HERE":
        return "(not set)"
    return f"{token[:8]}...{token[-4:]} (keep this secret!)"


print(f"Token: {_mask_token(IBM_TOKEN)}")


# ------------------------------------------------------------------
# IBM Quantum connection
# ------------------------------------------------------------------
def _try_connect(channel: str, token: str, instance=None):
    """Attempt a single QiskitRuntimeService connection and return the service."""
    from qiskit_ibm_runtime import QiskitRuntimeService  # type: ignore

    kwargs = dict(channel=channel, token=token)
    if instance is not None:
        kwargs["instance"] = instance
    return QiskitRuntimeService(**kwargs)


def _discover_instance(service) -> Optional[str]:
    """Return the name of the first available instance on *service*, or None."""
    try:
        instances = service.instances()
        if instances:
            return instances[0]
    except Exception:
        pass
    return None


def connect_ibm_service(token: str, preferred_channel: str, instance=None):
    """
    Connect to IBM Quantum, trying sensible channel / instance combinations.

    Resolution order
    ----------------
    1. preferred_channel + supplied instance (if any)
    2. preferred_channel + auto-discovered instance
    3. "ibm_quantum"  + auto-discovered instance
    4. "ibm_cloud"    + auto-discovered instance (no instance kwarg)

    Raises RuntimeError if every attempt fails.
    """
    # Normalise common mis-spellings / old names so the script is forgiving.
    _channel_aliases = {
        "ibm_quantum_platform": "ibm_quantum",
        "ibmq": "ibm_quantum",
        "quantum": "ibm_quantum",
        "cloud": "ibm_cloud",
    }
    preferred_channel = _channel_aliases.get(preferred_channel, preferred_channel)

    # Build an ordered list of (channel, instance) pairs to attempt.
    attempts: list = []  # list of (channel: str, instance: Optional[str])

    if instance:
        attempts.append((preferred_channel, instance))
    # Always also try without an explicit instance (SDK picks the default).
    attempts.append((preferred_channel, None))

    # If preferred channel is not "ibm_quantum" also try that as a fallback.
    if preferred_channel != "ibm_quantum":
        attempts.append(("ibm_quantum", None))
    # ibm_cloud as last resort.
    if preferred_channel != "ibm_cloud":
        attempts.append(("ibm_cloud", None))

    last_exc: Optional[Exception] = None
    for channel, inst in attempts:
        label = f"channel={channel}" + (f", instance={inst}" if inst else "")
        print(f"[IBM] Trying {label} ...")
        try:
            svc = _try_connect(channel, token, inst)
            # If no instance was given, try to discover one now.
            if inst is None:
                discovered = _discover_instance(svc)
                if discovered:
                    print(f"[IBM] Auto-discovered instance: {discovered}")
                    svc = _try_connect(channel, token, discovered)
            print(f"[IBM] Connected successfully (channel={channel})")
            return svc
        except Exception as exc:  # noqa: BLE001
            print(f"[IBM] '{channel}' failed: {exc}")
            last_exc = exc

    raise RuntimeError(f"Could not connect to IBM Quantum: {last_exc}")


# ------------------------------------------------------------------
# Hamiltonian
# ------------------------------------------------------------------
symbols = ["H", "H"]
coordinates = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.4])
H_h2, qubits = qml.qchem.molecular_hamiltonian(symbols, coordinates,
                                                charge=0, mult=1,
                                                basis='sto-3g')

if hasattr(H_h2, 'pauli_rep') and H_h2.pauli_rep is not None:
    num_terms = len(H_h2.pauli_rep)
elif hasattr(H_h2, 'ops'):
    num_terms = len(H_h2.ops)
elif hasattr(H_h2, 'operands'):
    num_terms = len(H_h2.operands)
else:
    num_terms = "N/A"

print(f"\n[VQE] H2 Hamiltonian: {qubits} qubits, {num_terms} Pauli terms")


# ------------------------------------------------------------------
# Ansatz
# ------------------------------------------------------------------
def hardware_efficient_ansatz(params, n_qubits: int):
    """Hardware-efficient ansatz: RY+RZ layers with CNOT entanglement."""
    n_layers = len(params)
    for layer in range(n_layers):
        for i in range(n_qubits):
            qml.RY(params[layer][i], wires=i)
            qml.RZ(params[layer][i + n_qubits], wires=i)
        for i in range(n_qubits - 1):
            qml.CNOT(wires=[i, i + 1])
        if n_qubits > 2:
            qml.CNOT(wires=[n_qubits - 1, 0])


# ------------------------------------------------------------------
# Device selection
# ------------------------------------------------------------------
def build_device(service=None):
    """
    Return a PennyLane device.

    If *service* is a connected QiskitRuntimeService, prefer an IBM backend.
    Otherwise fall back to 'default.qubit'.
    """
    if service is not None:
        try:
            backend_name = IBM_BACKEND or None
            if backend_name:
                backend = service.backend(backend_name)
            else:
                backend = service.least_busy(operational=True, simulator=True,
                                             min_num_qubits=qubits)
            print(f"[IBM] Using backend: {backend.name}")
            dev = qml.device('qiskit.ibmq', wires=qubits, backend=backend.name,
                             provider=service)
            return dev, backend.name
        except Exception as exc:  # noqa: BLE001
            print(f"[IBM] Could not attach PennyLane to IBM backend: {exc}")
            print("[IBM] Falling back to local simulation.")

    print("[Local] Using 'default.qubit' simulator.")
    return qml.device('default.qubit', wires=qubits), "default.qubit"


# ------------------------------------------------------------------
# VQE circuit
# ------------------------------------------------------------------
def make_vqe_circuit(dev):
    @qml.qnode(dev)
    def vqe_circuit(params, hamiltonian):
        # Hartree-Fock initial state |1100⟩ for H₂
        qml.PauliX(wires=0)
        qml.PauliX(wires=1)
        hardware_efficient_ansatz(params, qubits)
        return qml.expval(hamiltonian)
    return vqe_circuit


# ------------------------------------------------------------------
# Optimisation helpers
# ------------------------------------------------------------------
n_layers = 4
params_shape = (n_layers, 2 * qubits)


def run_vqe(vqe_circuit, init_params_flat):
    energy_history = []

    def callback(params):
        e = cost_function(params, vqe_circuit)
        energy_history.append(float(e))
        if len(energy_history) % 10 == 0:
            print(f"  Iteration {len(energy_history):3d}: Energy = {e:.8f} Ha")

    result = minimize(
        lambda p: cost_function(p, vqe_circuit),
        init_params_flat,
        method='L-BFGS-B',
        jac='3-point',
        callback=callback,
        options={'maxiter': 300, 'ftol': 1e-12},
    )
    return result, energy_history


def cost_function(params_flat, vqe_circuit):
    params = params_flat.reshape(params_shape)
    return float(vqe_circuit(params, H_h2))


# ------------------------------------------------------------------
# Plotting helpers
# ------------------------------------------------------------------
def plot_convergence(energy_history, hf_energy, exact_energy=-1.137283):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].plot(energy_history, 'b-', linewidth=2, label='VQE Energy')
    axes[0].axhline(y=exact_energy, color='r', linestyle='--', linewidth=2,
                    label='Exact Energy')
    axes[0].axhline(y=hf_energy, color='g', linestyle='--', linewidth=2,
                    label='Hartree-Fock')
    axes[0].set_xlabel('Optimization Iteration', fontsize=12)
    axes[0].set_ylabel('Energy (Hartree)', fontsize=12)
    axes[0].set_title('VQE Energy Convergence', fontsize=14, fontweight='bold')
    axes[0].legend(fontsize=10)
    axes[0].grid(True, alpha=0.3)

    errors = [abs(e - exact_energy) for e in energy_history]
    axes[1].semilogy(errors, 'b-', linewidth=2)
    axes[1].set_xlabel('Optimization Iteration', fontsize=12)
    axes[1].set_ylabel('Absolute Error (Hartree)', fontsize=12)
    axes[1].set_title('VQE Error Convergence (Log Scale)', fontsize=14,
                      fontweight='bold')
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('vqe_convergence.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("Convergence plot saved to vqe_convergence.png")


def plot_pes(bond_lengths, vqe_energies, hf_energies):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(bond_lengths * 0.529, vqe_energies, 'bo-', linewidth=2,
            markersize=8, label='VQE')
    ax.plot(bond_lengths * 0.529, hf_energies, 'gs--', linewidth=2,
            markersize=8, label='Hartree-Fock')
    ax.set_xlabel('H–H Bond Length (Ångström)', fontsize=12)
    ax.set_ylabel('Energy (Hartree)', fontsize=12)
    ax.set_title('H₂ Potential Energy Surface', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('h2_potential_energy_surface.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("PES plot saved to h2_potential_energy_surface.png")


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
def main():
    # 1. Connect (or fall back gracefully)
    service = None
    if IBM_TOKEN and IBM_TOKEN != "YOUR_IBM_TOKEN_HERE":
        try:
            service = connect_ibm_service(IBM_TOKEN, IBM_CHANNEL, IBM_INSTANCE)
        except RuntimeError as exc:
            print(f"\n[ERROR] {exc}")
            print("[ERROR] Check IBM_TOKEN and IBM_CHANNEL in the CONFIG section at the top.")
            if not FALLBACK_TO_LOCAL:
                raise
            print("[WARN]  Continuing with local simulation (FALLBACK_TO_LOCAL=True).\n")
    else:
        print("[WARN] IBM_TOKEN not configured – using local simulation.\n")

    # 2. Build device
    dev, backend_name = build_device(service)
    vqe_circuit = make_vqe_circuit(dev)

    # 3. Draw circuit
    np.random.seed(42)
    init_params = np.random.random(params_shape) * 0.1
    print("\nVQE Circuit Structure:")
    print(qml.draw(vqe_circuit)(init_params, H_h2))

    # 4. Initial energy
    init_energy = cost_function(init_params.flatten(), vqe_circuit)
    print(f"\nInitial energy: {init_energy:.8f} Hartree")
    print("\nStarting VQE optimisation...\n")

    # 5. Optimise
    result, energy_history = run_vqe(vqe_circuit, init_params.flatten())
    vqe_energy = result.fun

    # 6. HF reference
    mol = qml.qchem.Molecule(symbols,
                              np_std.array(coordinates.reshape(-1, 3)),
                              charge=0, mult=1, basis_name='sto-3g')
    hf_energy = float(qml.qchem.scf(mol)()[0][0])

    exact_energy = -1.137283

    print(f"\n{'='*50}")
    print("VQE Optimisation Complete!")
    print(f"Final Ground State Energy : {vqe_energy:.8f} Hartree")
    print(f"{'='*50}")

    vqe_error = abs(vqe_energy - exact_energy)
    hf_error  = abs(hf_energy  - exact_energy)

    print("\n" + "="*60)
    print("COMPARISON WITH CLASSICAL METHODS")
    print("="*60)
    print(f"\nExact Ground State Energy (FCI): {exact_energy:.8f} Ha")
    print(f"Hartree-Fock Energy (Classical):  {hf_energy:.8f} Ha")
    print(f"VQE Energy (Quantum):             {vqe_energy:.8f} Ha")
    print(f"\nAbsolute Errors:")
    print(f"  Hartree-Fock Error: {hf_error:.8f} Ha  ({hf_error*627.5:.4f} kcal/mol)")
    print(f"  VQE Error:          {vqe_error:.8f} Ha  ({vqe_error*627.5:.4f} kcal/mol)")
    if hf_error > 0:
        print(f"\nImprovement over HF: {((hf_error - vqe_error)/hf_error * 100):.2f}%")

    if vqe_error < hf_error:
        print("\n✓ VQE achieved better accuracy than Hartree-Fock!")
    else:
        print("\n⚠ VQE needs more optimisation layers or iterations.")

    # 7. Convergence plot
    plot_convergence(energy_history, hf_energy, exact_energy)

    # 8. Potential Energy Surface scan
    bond_lengths = np.linspace(0.5, 4.5, 15)
    vqe_energies, hf_energies_pes = [], []

    print("\nScanning potential energy surface...\n")
    current_params = result.x.copy()

    for bond_length in bond_lengths:
        coords = np.array([0.0, 0.0, 0.0, 0.0, 0.0, bond_length])
        H_scan, _ = qml.qchem.molecular_hamiltonian(symbols, coords, basis='sto-3g')

        @qml.qnode(dev)
        def _circuit_scan(p, H=H_scan):
            qml.PauliX(wires=0)
            qml.PauliX(wires=1)
            hardware_efficient_ansatz(p, qubits)
            return qml.expval(H)

        res_scan = minimize(
            lambda p: float(_circuit_scan(p.reshape(params_shape))),
            current_params,
            method='L-BFGS-B',
            jac='3-point',
            options={'maxiter': 200, 'ftol': 1e-12},
        )
        energy_vqe = res_scan.fun
        current_params = res_scan.x.copy()
        vqe_energies.append(energy_vqe)

        mol_scan = qml.qchem.Molecule(symbols,
                                       np_std.array(coords.reshape(-1, 3)),
                                       basis_name='sto-3g')
        energy_hf = float(qml.qchem.scf(mol_scan)()[0][0])
        hf_energies_pes.append(energy_hf)

        print(f"  R = {bond_length:.2f} Bohr ({bond_length*0.529:.2f} Å) | "
              f"VQE: {energy_vqe:.6f} Ha | HF: {energy_hf:.6f} Ha")

    plot_pes(bond_lengths, vqe_energies, hf_energies_pes)

    opt_idx = int(np.argmin(vqe_energies))
    print(f"\nOptimal bond length : {bond_lengths[opt_idx] * 0.529:.3f} Å")
    print(f"Minimum VQE energy  : {min(vqe_energies):.8f} Ha")


if __name__ == "__main__":
    main()
