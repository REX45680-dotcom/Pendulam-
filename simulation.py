"""
Computational Analysis for Research Paper:
"Deterministic Chaos in the Double Pendulum: Lyapunov Exponent Analysis 
 and Phase Space Topology"

Author: [Your Name]
Description: This script performs numerical simulations of the double pendulum
system, computes Lyapunov exponents, generates phase space portraits, 
Poincaré sections, and bifurcation diagrams.

Requirements: numpy, scipy, matplotlib
Install: pip install numpy scipy matplotlib
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os

# ============================================================================
# PHYSICAL CONSTANTS AND SYSTEM PARAMETERS
# ============================================================================

class DoublePendulumParams:
    """Parameters for the double pendulum system."""
    def __init__(self, m1=1.0, m2=1.0, L1=1.0, L2=1.0, g=9.81):
        self.m1 = m1  # Mass of pendulum 1 (kg)
        self.m2 = m2  # Mass of pendulum 2 (kg)
        self.L1 = L1  # Length of pendulum 1 (m)
        self.L2 = L2  # Length of pendulum 2 (m)
        self.g = g    # Gravitational acceleration (m/s^2)


# ============================================================================
# EQUATIONS OF MOTION
# ============================================================================

def double_pendulum_equations(t, state, params):
    """
    Full nonlinear equations of motion for the double pendulum.
    
    Uses Lagrangian mechanics derivation:
    The Lagrangian L = T - V where T is kinetic energy and V is potential energy.
    
    State vector: [theta1, omega1, theta2, omega2]
    where omega_i = d(theta_i)/dt
    
    Returns: [dtheta1/dt, domega1/dt, dtheta2/dt, domega2/dt]
    """
    theta1, omega1, theta2, omega2 = state
    m1, m2, L1, L2, g = params.m1, params.m2, params.L1, params.L2, params.g
    
    delta = theta2 - theta1
    
    # Denominators from the mass matrix inversion
    den1 = (m1 + m2) * L1 - m2 * L1 * np.cos(delta)**2
    den2 = (L2 / L1) * den1
    
    # Angular accelerations derived from Euler-Lagrange equations
    dtheta1_dt = omega1
    dtheta2_dt = omega2
    
    domega1_dt = (m2 * L1 * omega1**2 * np.sin(delta) * np.cos(delta) +
                  m2 * g * np.sin(theta2) * np.cos(delta) +
                  m2 * L2 * omega2**2 * np.sin(delta) -
                  (m1 + m2) * g * np.sin(theta1)) / den1
    
    domega2_dt = (-m2 * L2 * omega2**2 * np.sin(delta) * np.cos(delta) +
                  (m1 + m2) * g * np.sin(theta1) * np.cos(delta) -
                  (m1 + m2) * L1 * omega1**2 * np.sin(delta) -
                  (m1 + m2) * g * np.sin(theta2)) / den2
    
    return [dtheta1_dt, domega1_dt, dtheta2_dt, domega2_dt]


def compute_energy(state, params):
    """
    Compute total mechanical energy E = T + V of the double pendulum.
    Used to verify numerical accuracy (energy should be conserved).
    """
    theta1, omega1, theta2, omega2 = state.T if state.ndim > 1 else state
    m1, m2, L1, L2, g = params.m1, params.m2, params.L1, params.L2, params.g
    
    # Kinetic energy
    T = (0.5 * (m1 + m2) * L1**2 * omega1**2 +
         0.5 * m2 * L2**2 * omega2**2 +
         m2 * L1 * L2 * omega1 * omega2 * np.cos(theta1 - theta2))
    
    # Potential energy (reference: both pendulums hanging straight down)
    V = (-(m1 + m2) * g * L1 * np.cos(theta1) -
          m2 * g * L2 * np.cos(theta2))
    
    return T + V


# ============================================================================
# SIMULATION ENGINE
# ============================================================================

def simulate_double_pendulum(params, initial_state, t_span, dt=0.001):
    """
    Simulate double pendulum using RK45 (Runge-Kutta 4th/5th order).
    
    Args:
        params: DoublePendulumParams instance
        initial_state: [theta1_0, omega1_0, theta2_0, omega2_0]
        t_span: (t_start, t_end)
        dt: Maximum time step
    
    Returns:
        t: Time array
        solution: State array [theta1(t), omega1(t), theta2(t), omega2(t)]
    """
    t_eval = np.arange(t_span[0], t_span[1], dt)
    
    sol = solve_ivp(
        double_pendulum_equations,
        t_span,
        initial_state,
        args=(params,),
        method='RK45',
        t_eval=t_eval,
        rtol=1e-10,
        atol=1e-12,
        max_step=dt
    )
    
    return sol.t, sol.y


# ============================================================================
# LYAPUNOV EXPONENT COMPUTATION
# ============================================================================

def compute_lyapunov_exponents(params, initial_state, t_total=100.0, dt=0.01,
                                renorm_steps=100):
    """
    Compute the Maximal Lyapunov Exponent (MLE) using the tangent space method.
    
    The MLE quantifies the rate of exponential divergence of nearby trajectories:
        λ = lim(t→∞) (1/t) * ln(|δx(t)| / |δx(0)|)
    
    λ > 0 indicates chaos
    λ ≈ 0 indicates quasi-periodic motion
    λ < 0 indicates stable fixed point
    
    Method: 
    1. Evolve a reference trajectory and a nearby perturbed trajectory
    2. Periodically renormalize the perturbation to prevent overflow
    3. Accumulate the logarithm of the stretching factors
    """
    epsilon = 1e-8  # Initial perturbation magnitude
    
    # Reference trajectory
    ref_state = np.array(initial_state, dtype=float)
    
    # Perturbed trajectory (small perturbation in theta1)
    pert_state = ref_state.copy()
    pert_state[0] += epsilon
    
    lyap_sum = 0.0
    lyap_history = []
    t_history = []
    n_steps = int(t_total / (dt * renorm_steps))
    
    for i in range(1, n_steps + 1):
        t_span_step = (0, dt * renorm_steps)
        t_eval_step = np.array([dt * renorm_steps])
        
        # Evolve reference
        sol_ref = solve_ivp(
            double_pendulum_equations, t_span_step, ref_state,
            args=(params,), method='RK45', t_eval=t_eval_step,
            rtol=1e-10, atol=1e-12
        )
        ref_state = sol_ref.y[:, -1]
        
        # Evolve perturbed
        sol_pert = solve_ivp(
            double_pendulum_equations, t_span_step, pert_state,
            args=(params,), method='RK45', t_eval=t_eval_step,
            rtol=1e-10, atol=1e-12
        )
        pert_state = sol_pert.y[:, -1]
        
        # Compute separation
        delta = pert_state - ref_state
        dist = np.linalg.norm(delta)
        
        if dist > 0:
            lyap_sum += np.log(dist / epsilon)
        
        # Renormalize perturbation
        pert_state = ref_state + epsilon * (delta / dist)
        
        current_time = i * dt * renorm_steps
        lyap_history.append(lyap_sum / current_time)
        t_history.append(current_time)
    
    mle = lyap_sum / t_total
    return mle, np.array(t_history), np.array(lyap_history)


# ============================================================================
# SENSITIVITY TO INITIAL CONDITIONS
# ============================================================================

def sensitivity_analysis(params, base_state, perturbations, t_span, dt=0.001):
    """
    Demonstrate sensitive dependence on initial conditions — 
    the hallmark of deterministic chaos.
    
    Simulates multiple trajectories with slightly different initial conditions
    and tracks their divergence over time.
    """
    trajectories = []
    
    for eps in perturbations:
        perturbed = base_state.copy()
        perturbed[0] += eps  # Perturb theta1
        t, y = simulate_double_pendulum(params, perturbed, t_span, dt)
        trajectories.append((t, y, eps))
    
    return trajectories


# ============================================================================
# POINCARÉ SECTION
# ============================================================================

def compute_poincare_section(t, state, section_var_idx=0, section_value=0.0,
                              velocity_sign=1):
    """
    Compute Poincaré section by finding intersections of the trajectory
    with the hyperplane θ_i = section_value (with specified velocity sign).
    
    This reduces the continuous flow to a discrete map, revealing the 
    geometric structure of the attractor.
    """
    var = state[section_var_idx]
    vel = state[section_var_idx + 1]
    
    crossings_x = []
    crossings_y = []
    
    for i in range(len(var) - 1):
        # Detect zero crossing with correct velocity direction
        if ((var[i] - section_value) * (var[i+1] - section_value) < 0 and
            velocity_sign * vel[i] > 0):
            # Linear interpolation for better accuracy
            frac = (section_value - var[i]) / (var[i+1] - var[i])
            
            # Record the other pair of phase space variables
            other_var_idx = 2 if section_var_idx == 0 else 0
            x_cross = state[other_var_idx][i] + frac * (state[other_var_idx][i+1] - state[other_var_idx][i])
            y_cross = state[other_var_idx + 1][i] + frac * (state[other_var_idx + 1][i+1] - state[other_var_idx + 1][i])
            
            crossings_x.append(x_cross)
            crossings_y.append(y_cross)
    
    return np.array(crossings_x), np.array(crossings_y)


# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

def plot_phase_space_portraits(t, state, params, save_path):
    """Generate phase space portraits for both pendulums."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Phase Space Portraits of the Double Pendulum', 
                 fontsize=14, fontweight='bold')
    
    # Pendulum 1
    axes[0].plot(state[0], state[1], linewidth=0.3, color='#2196F3', alpha=0.7)
    axes[0].set_xlabel(r'$\theta_1$ (rad)', fontsize=12)
    axes[0].set_ylabel(r'$\dot{\theta}_1$ (rad/s)', fontsize=12)
    axes[0].set_title('Pendulum 1: Phase Portrait', fontsize=12)
    axes[0].grid(True, alpha=0.3)
    
    # Pendulum 2
    axes[1].plot(state[2], state[3], linewidth=0.3, color='#FF5722', alpha=0.7)
    axes[1].set_xlabel(r'$\theta_2$ (rad)', fontsize=12)
    axes[1].set_ylabel(r'$\dot{\theta}_2$ (rad/s)', fontsize=12)
    axes[1].set_title('Pendulum 2: Phase Portrait', fontsize=12)
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Phase space portraits saved: {save_path}")


def plot_trajectory_comparison(t, state, params, save_path):
    """Plot angular displacement time series."""
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    fig.suptitle('Angular Displacement Time Series', fontsize=14, fontweight='bold')
    
    axes[0].plot(t, np.degrees(state[0]), linewidth=0.5, color='#2196F3')
    axes[0].set_ylabel(r'$\theta_1$ (degrees)', fontsize=12)
    axes[0].set_title('Pendulum 1', fontsize=12)
    axes[0].grid(True, alpha=0.3)
    
    axes[1].plot(t, np.degrees(state[2]), linewidth=0.5, color='#FF5722')
    axes[1].set_ylabel(r'$\theta_2$ (degrees)', fontsize=12)
    axes[1].set_xlabel('Time (s)', fontsize=12)
    axes[1].set_title('Pendulum 2', fontsize=12)
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Trajectory comparison saved: {save_path}")


def plot_sensitivity(trajectories, save_path):
    """Visualize sensitive dependence on initial conditions."""
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))
    fig.suptitle('Sensitive Dependence on Initial Conditions', 
                 fontsize=14, fontweight='bold')
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(trajectories)))
    
    for (t, y, eps), color in zip(trajectories, colors):
        label = f'Δθ₁ = {eps:.1e} rad'
        axes[0].plot(t, np.degrees(y[0]), linewidth=0.5, color=color, label=label)
        axes[1].plot(t, np.degrees(y[2]), linewidth=0.5, color=color, label=label)
    
    axes[0].set_ylabel(r'$\theta_1$ (degrees)', fontsize=12)
    axes[0].set_title('Pendulum 1 — Diverging Trajectories', fontsize=12)
    axes[0].legend(fontsize=8, loc='upper right')
    axes[0].grid(True, alpha=0.3)
    
    axes[1].set_ylabel(r'$\theta_2$ (degrees)', fontsize=12)
    axes[1].set_xlabel('Time (s)', fontsize=12)
    axes[1].set_title('Pendulum 2 — Diverging Trajectories', fontsize=12)
    axes[1].legend(fontsize=8, loc='upper right')
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Sensitivity analysis saved: {save_path}")


def plot_lyapunov_convergence(t_hist, lyap_hist, mle, save_path):
    """Plot convergence of the Maximal Lyapunov Exponent."""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.plot(t_hist, lyap_hist, linewidth=1.5, color='#9C27B0')
    ax.axhline(y=mle, color='red', linestyle='--', linewidth=1, 
               label=f'MLE = {mle:.4f} bits/s')
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)
    
    ax.set_xlabel('Time (s)', fontsize=12)
    ax.set_ylabel('Lyapunov Exponent (bits/s)', fontsize=12)
    ax.set_title('Convergence of the Maximal Lyapunov Exponent', 
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    
    # Add chaos/order annotation
    if mle > 0:
        ax.annotate(f'λ > 0 → CHAOTIC REGIME', xy=(0.7, 0.9),
                   xycoords='axes fraction', fontsize=12, color='red',
                   fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow'))
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Lyapunov convergence saved: {save_path}")


def plot_poincare_section(crossings_x, crossings_y, save_path):
    """Plot Poincaré section."""
    fig, ax = plt.subplots(figsize=(8, 8))
    
    ax.scatter(crossings_x, crossings_y, s=0.5, color='#E91E63', alpha=0.6)
    ax.set_xlabel(r'$\theta_2$ (rad)', fontsize=12)
    ax.set_ylabel(r'$\dot{\theta}_2$ (rad/s)', fontsize=12)
    ax.set_title('Poincaré Section (θ₁ = 0, dθ₁/dt > 0)', 
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('auto')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Poincaré section saved: {save_path}")


def plot_energy_conservation(t, state, params, save_path):
    """Plot total energy over time to verify numerical accuracy."""
    energies = []
    for i in range(len(t)):
        s = state[:, i]
        E = compute_energy(s, params)
        energies.append(E)
    energies = np.array(energies)
    
    E0 = energies[0]
    relative_error = np.abs((energies - E0) / E0) if E0 != 0 else np.abs(energies - E0)
    
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))
    fig.suptitle('Energy Conservation Verification', fontsize=14, fontweight='bold')
    
    axes[0].plot(t, energies, linewidth=0.5, color='#4CAF50')
    axes[0].set_ylabel('Total Energy E (J)', fontsize=12)
    axes[0].set_title('Total Mechanical Energy vs Time', fontsize=12)
    axes[0].grid(True, alpha=0.3)
    
    axes[1].semilogy(t, relative_error, linewidth=0.5, color='#FF9800')
    axes[1].set_ylabel('Relative Energy Error |ΔE/E₀|', fontsize=12)
    axes[1].set_xlabel('Time (s)', fontsize=12)
    axes[1].set_title('Numerical Integration Accuracy', fontsize=12)
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Energy conservation saved: {save_path}")


def plot_lyapunov_map(params, theta1_range, theta2_range, resolution, save_path):
    """
    Generate a Lyapunov exponent map over initial condition space.
    This reveals the fractal boundary between regular and chaotic regions.
    """
    theta1_vals = np.linspace(*theta1_range, resolution)
    theta2_vals = np.linspace(*theta2_range, resolution)
    lyap_map = np.zeros((resolution, resolution))
    
    total = resolution * resolution
    print(f"  Computing Lyapunov map ({resolution}x{resolution} = {total} points)...")
    
    for i, th1 in enumerate(theta1_vals):
        if i % 5 == 0:
            print(f"    Row {i+1}/{resolution}...")
        for j, th2 in enumerate(theta2_vals):
            initial = [th1, 0.0, th2, 0.0]
            mle, _, _ = compute_lyapunov_exponents(
                params, initial, t_total=20.0, dt=0.01, renorm_steps=50
            )
            lyap_map[j, i] = mle
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    im = ax.imshow(lyap_map, extent=[*theta1_range, *theta2_range],
                   origin='lower', cmap='inferno', aspect='auto',
                   vmin=-0.5, vmax=2.0)
    
    cbar = plt.colorbar(im, ax=ax, label='Maximal Lyapunov Exponent (bits/s)')
    ax.set_xlabel(r'$\theta_1$ (rad)', fontsize=12)
    ax.set_ylabel(r'$\theta_2$ (rad)', fontsize=12)
    ax.set_title('Lyapunov Exponent Map: Order vs Chaos in Initial Condition Space',
                 fontsize=13, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Lyapunov map saved: {save_path}")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Run all analyses and generate figures for the research paper."""
    
    output_dir = os.path.dirname(os.path.abspath(__file__))
    figures_dir = os.path.join(output_dir, "figures")
    os.makedirs(figures_dir, exist_ok=True)
    
    # System parameters (equal mass, equal length)
    params = DoublePendulumParams(m1=1.0, m2=1.0, L1=1.0, L2=1.0, g=9.81)
    
    # ---- Chaotic Initial Conditions ----
    # Large initial angles → chaotic regime
    theta1_0 = np.radians(120)  # 120 degrees
    theta2_0 = np.radians(150)  # 150 degrees
    omega1_0 = 0.0
    omega2_0 = 0.0
    initial_chaotic = [theta1_0, omega1_0, theta2_0, omega2_0]
    
    print("=" * 70)
    print("  DOUBLE PENDULUM CHAOS ANALYSIS")
    print("  Research Paper Computational Supplement")
    print("=" * 70)
    
    # ---- 1. Main Simulation ----
    print("\n[1/6] Running main simulation (t = 0 to 60s)...")
    t_span = (0, 60)
    t, state = simulate_double_pendulum(params, initial_chaotic, t_span, dt=0.001)
    print(f"  ✓ Simulation complete: {len(t)} time steps")
    
    # ---- 2. Phase Space Portraits ----
    print("\n[2/6] Generating phase space portraits...")
    plot_phase_space_portraits(t, state, params, 
                               os.path.join(figures_dir, "fig1_phase_space.png"))
    plot_trajectory_comparison(t, state, params,
                               os.path.join(figures_dir, "fig2_time_series.png"))
    
    # ---- 3. Sensitivity Analysis ----
    print("\n[3/6] Running sensitivity analysis...")
    perturbations = [0, 1e-10, 1e-8, 1e-6, 1e-4]
    trajs = sensitivity_analysis(params, initial_chaotic, perturbations, 
                                  (0, 30), dt=0.001)
    plot_sensitivity(trajs, os.path.join(figures_dir, "fig3_sensitivity.png"))
    
    # ---- 4. Lyapunov Exponent ----
    print("\n[4/6] Computing Maximal Lyapunov Exponent...")
    mle, t_hist, lyap_hist = compute_lyapunov_exponents(
        params, initial_chaotic, t_total=200.0, dt=0.01, renorm_steps=100
    )
    print(f"  ✓ Maximal Lyapunov Exponent (MLE) = {mle:.6f} bits/s")
    if mle > 0:
        print(f"  → CHAOTIC REGIME CONFIRMED (λ > 0)")
        print(f"  → Predictability horizon ≈ {1/mle:.2f} seconds")
    plot_lyapunov_convergence(t_hist, lyap_hist, mle,
                              os.path.join(figures_dir, "fig4_lyapunov.png"))
    
    # ---- 5. Poincaré Section ----
    print("\n[5/6] Computing Poincaré section...")
    t_long, state_long = simulate_double_pendulum(params, initial_chaotic, 
                                                    (0, 500), dt=0.001)
    cx, cy = compute_poincare_section(t_long, state_long, 
                                       section_var_idx=0, section_value=0.0,
                                       velocity_sign=1)
    print(f"  ✓ Found {len(cx)} crossings")
    plot_poincare_section(cx, cy, os.path.join(figures_dir, "fig5_poincare.png"))
    
    # ---- 6. Energy Conservation Check ----
    print("\n[6/6] Verifying energy conservation...")
    plot_energy_conservation(t, state, params,
                             os.path.join(figures_dir, "fig6_energy.png"))
    
    # ---- Summary Statistics ----
    E_initial = compute_energy(state[:, 0], params)
    E_final = compute_energy(state[:, -1], params)
    
    print("\n" + "=" * 70)
    print("  RESULTS SUMMARY")
    print("=" * 70)
    print(f"  Initial conditions: θ₁={np.degrees(theta1_0):.0f}°, θ₂={np.degrees(theta2_0):.0f}°")
    print(f"  Maximal Lyapunov Exponent: λ = {mle:.6f} bits/s")
    print(f"  Predictability horizon: τ ≈ {1/mle:.2f} s")
    print(f"  Initial energy: E₀ = {E_initial:.6f} J")
    print(f"  Final energy:   E_f = {E_final:.6f} J")
    print(f"  Relative energy error: |ΔE/E₀| = {abs((E_final-E_initial)/E_initial):.2e}")
    print(f"  Poincaré crossings: {len(cx)}")
    print(f"  Figures saved to: {figures_dir}")
    print("=" * 70)
    
    # ---- Optional: Lyapunov Map (computationally intensive) ----
    print("\n[OPTIONAL] Lyapunov exponent map (30x30 grid)...")
    print("  (This may take a few minutes...)")
    plot_lyapunov_map(params, 
                      theta1_range=(0.1, np.pi),
                      theta2_range=(0.1, np.pi),
                      resolution=30,
                      save_path=os.path.join(figures_dir, "fig7_lyapunov_map.png"))
    
    print("\n✅ All analyses complete! Figures ready for the research paper.")


if __name__ == "__main__":
    main()
