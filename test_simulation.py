"""
Automated Verification & Unit Tests for Double Pendulum Simulation
Matches Terminal-Bench / Benchmark scorer standards.
"""

import pytest
import numpy as np
from simulation import (
    DoublePendulumParams,
    double_pendulum_equations,
    compute_energy,
    simulate_double_pendulum,
    compute_lyapunov_exponents
)

@pytest.fixture
def standard_params():
    return DoublePendulumParams(m1=1.0, m2=1.0, L1=1.0, L2=1.0, g=9.81)

def test_equilibrium_derivative_zero(standard_params):
    """The hanging equilibrium (0, 0, 0, 0) must have zero time derivatives."""
    state = [0.0, 0.0, 0.0, 0.0]
    derivs = double_pendulum_equations(0.0, state, standard_params)
    assert np.allclose(derivs, [0.0, 0.0, 0.0, 0.0], atol=1e-12)

def test_energy_conservation(standard_params):
    """Energy should be conserved within numerical tolerance over short integration."""
    initial_state = [np.radians(45), 0.0, np.radians(30), 0.0]
    t, state = simulate_double_pendulum(standard_params, initial_state, (0, 5), dt=0.002)
    
    e_start = compute_energy(state[:, 0], standard_params)
    e_end = compute_energy(state[:, -1], standard_params)
    rel_error = abs((e_end - e_start) / e_start)
    
    assert rel_error < 1e-6, f"Energy drift too high: {rel_error}"

def test_lyapunov_positive_for_large_angle(standard_params):
    """Large amplitude motion must produce positive Lyapunov exponent (chaos)."""
    initial_state = [np.radians(120), 0.0, np.radians(150), 0.0]
    mle, _, _ = compute_lyapunov_exponents(
        standard_params, initial_state, t_total=15.0, dt=0.01, renorm_steps=50
    )
    assert mle > 0.5, f"Expected chaotic regime (mle > 0.5), got {mle}"
