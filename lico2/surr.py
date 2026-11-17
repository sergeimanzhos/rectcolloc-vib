"""Degree-4 polynomial surrogate of the DFT points (used by run_rc.py for the centre selection V <= Ec).
Reconstructed 14 Nov 2026: surr.py was not in the repository; run_rc.py calls polyfit(q, E) -> (callable, rms cm-1).
Code copied from dft2.polyfit (all monomials of total degree <= 4 in the six normal coordinates, least squares)."""
import itertools, numpy as np
import coords as c
H = c.H2CM
EXPS = np.array([e for d in range(5) for e in itertools.product(range(d + 1), repeat=c.NQ) if sum(e) == d])
def Phi(Q):
    Q = np.atleast_2d(Q); out = np.ones((len(Q), len(EXPS)))
    for k in range(c.NQ):
        pw = np.stack([Q[:, k] ** n for n in range(5)], 1); out *= pw[:, EXPS[:, k]]
    return out
def polyfit(q, E):
    coef = np.linalg.lstsq(Phi(q), E, rcond=None)[0]
    rms = np.sqrt(np.mean((Phi(q) @ coef - E) ** 2)) * H
    return (lambda Q: Phi(Q) @ coef), rms
