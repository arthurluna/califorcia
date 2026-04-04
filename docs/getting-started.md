# Getting Started

This guide shows the typical workflow for defining a Casimir-Lifshitz system and computing observables with `califorcia`.

## 1. Import The Main API

```python
from califorcia import system
from califorcia.materials import gold, vacuum
```

The package exposes one main public class: `system`.

## 2. Define The Geometry

Create a system by specifying:

- `T`: temperature in kelvin
- `d`: plate separation in meters
- `matL`: material of the left plane
- `matR`: material of the right plane
- `matm`: material of the intervening medium

Example:

```python
T = 300.0
d = 1.0e-6

s = system(T, d, gold, gold, vacuum)
```

For homogeneous half-spaces, `matL` and `matR` can be a single material object or module.

## 3. Compute Observables

The main convenience methods are:

```python
energy = s.energy()
pressure = s.pressure()
gradient = s.pressuregradient()
```

These correspond to:

- free energy per unit area
- pressure
- pressure gradient

All results are returned in SI units.

## 4. Choose The Frequency-Summation Method

At finite temperature, the package supports two summation methods:

- `fs="psd"`: Pade spectrum decomposition
- `fs="msd"`: direct Matsubara summation

Example:

```python
pressure = s.pressure(fs="psd", epsrel=1e-8)
```

Parameters:

- `epsrel`: target relative precision for the frequency summation
- `N`: optional number of summation terms

If `N` is omitted, the package chooses it automatically.

## 5. Zero-Temperature And High-Temperature Cases

If `T == 0`, the package performs a continuous frequency integration.

At finite temperature, you can request only the zero-frequency contribution through:

```python
s.energy(ht_limit=True)
```

This is useful when studying the high-temperature limit.

## 6. Add Coatings

Coated planes are represented by lists. The first material is the layer facing the medium, and the last material is the substrate half-space.

```python
from califorcia.materials import gold, teflon, vacuum

s = system(
    300.0,
    1.0e-6,
    [teflon, gold],
    gold,
    vacuum,
    deltaL=[50e-9],
)
```

Rules:

- `len(matL) == len(deltaL) + 1`
- `len(matR) == len(deltaR) + 1`
- up to three coating layers are supported on each side

## 7. Define A Custom Material

A custom material must provide:

- `materialclass`
- `epsilon(xi)`

Example:

```python
class UserMaterial:
    def __init__(self):
        self.materialclass = "dielectric"

    def epsilon(self, xi):
        wj = 1.911e15
        cj = 1.282
        return 1.0 + cj * wj**2 / (wj**2 + xi**2)
```

Then use it like any built-in material:

```python
u = UserMaterial()
s = system(300.0, 100e-9, gold, u, vacuum)
```

See the materials guide for the supported `materialclass` values.
