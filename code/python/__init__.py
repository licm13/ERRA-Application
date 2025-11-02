"""ERRA Python工具集。

This package provides a lightweight Python implementation of the
Ensemble Rainfall-Runoff Analysis (ERRA) workflow that accompanies the
original R scripts published by James Kirchner.  The Python port focuses
on the core regression engine for estimating Runoff Response
Distributions (RRDs) and related diagnostics, enabling integration into
modern Python-based hydrologic analyses.
"""

from .erra import ERRAResult, erra

__all__ = ["ERRAResult", "erra"]
