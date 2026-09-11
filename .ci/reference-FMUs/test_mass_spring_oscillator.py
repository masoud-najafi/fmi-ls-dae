"""Simulation tests for the MassSpringOscillator reference FMU.

Run with:
    MASS_SPRING_OSCILLATOR_FMU=path/to/MassSpringOscillator.fmu pytest test_mass_spring_oscillator.py -v
"""

import os

import numpy as np
import pytest
from fmpy import simulate_fmu

FMU = os.environ.get('MASS_SPRING_OSCILLATOR_FMU', 'MassSpringOscillator.fmu')

STOP_TIME = 20.0


def test_ode_mode():
    """FMU reproduces the analytical solution in ODE (Model Exchange) mode.

    A non-DAE-aware importer resolves the spring force analytically to F = -k*x
    and integrates the remaining ODE.  With the default parameters m = k = 1 and
    the initial state (x, v) = (1, 0) the analytical solution is

        x(t) =  cos(t)
        v(t) = -sin(t)
    """
    result = simulate_fmu(
        FMU,
        start_time=0.0,
        stop_time=STOP_TIME,
        output_interval=0.1,
        solver='CVode',
        relative_tolerance=1e-8,
        output=['x', 'v', 'E', '__residual0'],
    )

    assert result['time'][-1] >= STOP_TIME - 1e-6, \
        f"simulation stopped early at t={result['time'][-1]}"

    t = result['time']
    atol = 1e-4

    assert np.allclose(result['x'], np.cos(t), atol=atol), \
        f"x: max deviation from cos(t) is {np.max(np.abs(result['x'] - np.cos(t)))}"
    assert np.allclose(result['v'], -np.sin(t), atol=atol), \
        f"v: max deviation from -sin(t) is {np.max(np.abs(result['v'] + np.sin(t)))}"

    # F = -k*x holds exactly in ODE mode, so the residual is identically zero
    assert np.all(result['__residual0'] == 0.0), \
        f"__residual0 is not zero, max |residual| = {np.max(np.abs(result['__residual0']))}"


def test_invariant_is_preserved_by_an_accurate_solver():
    """The invariant E stays at its initial value when the integration is accurate."""
    result = simulate_fmu(
        FMU,
        start_time=0.0,
        stop_time=STOP_TIME,
        output_interval=0.1,
        solver='CVode',
        relative_tolerance=1e-8,
        output=['E'],
    )

    E = result['E']
    E0 = 0.5  # m/2*v0^2 + k/2*x0^2 with m = k = 1, x0 = 1, v0 = 0

    assert np.isclose(E[0], E0, atol=1e-12), f"E(0): expected {E0}, got {E[0]}"
    assert np.allclose(E, E0, atol=1e-4), \
        f"E drifted by {np.max(np.abs(E - E0))} despite a tight solver tolerance"


@pytest.mark.parametrize('h', [0.2, 0.1, 0.05])
def test_invariant_drifts_with_explicit_euler(h):
    """The invariant drifts exactly as predicted when integrated with explicit Euler.

    Substituting the explicit Euler step into E gives, for m = k = 1,

        E_{n+1} = E_n * (1 + h^2)

    so after N steps the energy has grown to E_0 * (1 + h^2)^N.  This is the
    numerical drift that motivates exposing invariants; see the "Invariants and
    Numerical Drift" appendix of the layered standard.
    """
    result = simulate_fmu(
        FMU,
        start_time=0.0,
        stop_time=STOP_TIME,
        output_interval=h,
        solver='Euler',
        step_size=h,
        output=['E'],
    )

    E = result['E']
    n_steps = int(round(STOP_TIME / h))
    expected = 0.5 * (1.0 + h * h) ** n_steps

    assert E[-1] > E[0], "explicit Euler is expected to increase the energy"
    assert np.isclose(E[-1], expected, rtol=1e-9), \
        f"E({STOP_TIME}) with h={h}: expected {expected}, got {E[-1]}"


def test_parameters():
    """Parameters m and k and the initial state can be changed by the importer."""
    m, k, x0 = 2.0, 8.0, 0.5
    omega = np.sqrt(k / m)

    result = simulate_fmu(
        FMU,
        start_time=0.0,
        stop_time=10.0,
        output_interval=0.1,
        solver='CVode',
        relative_tolerance=1e-8,
        start_values={'m': m, 'k': k, 'x': x0, 'v': 0.0},
        output=['x', 'E'],
    )

    t = result['time']

    assert np.allclose(result['x'], x0 * np.cos(omega * t), atol=1e-4), \
        f"x: max deviation is {np.max(np.abs(result['x'] - x0 * np.cos(omega * t)))}"
    assert np.isclose(result['E'][0], 0.5 * k * x0 ** 2, atol=1e-12), \
        f"E(0): expected {0.5 * k * x0 ** 2}, got {result['E'][0]}"


@pytest.mark.skip(reason="DAE-aware importer not yet available")
def test_dae_mode():
    """Simulate with a DAE-aware importer and verify the constraint and the invariant.

    A DAE-aware importer must:
    - read AlgebraicVariables and ModelStructure from
      extra/org.fmi-standard.fmi-ls-dae/fmi-ls-manifest.xml
    - solve 0 = F + k*x for the algebraic variable F at each step via fmi3SetFloat64
    - read the <Invariant> element and monitor E to detect drift

    TODO: replace the pass below with a real DAE importer call once one is
    available and assert, e.g.:
        assert np.all(np.abs(result['__residual0']) < 1e-8)
        assert np.all(np.abs(result['E'] - result['E'][0]) < 1e-6)
    """
    pass
