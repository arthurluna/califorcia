# Materials Guide

`califorcia` models materials through a minimal interface that supplies the dielectric response on the imaginary-frequency axis.

## Built-In Materials

The package includes predefined materials under `califorcia.materials`, including:

- `vacuum`
- `pec`
- `gold`
- `gold_drude`
- `gold_plasma`
- `aluminium`
- `teflon`
- `ethanol`
- `silica`
- `fused_silica`
- `aSiO2`
- `aAl2O3`
- `Si`
- `pdoped_Si`
- `SiC`
- `polystyrene`
- `sodalime`

Import them as:

```python
from califorcia.materials import gold, vacuum, teflon
```

## Material Interface

A material must expose:

- `materialclass`
- `epsilon(xi)`

Example:

```python
class UserMaterial:
    def __init__(self):
        self.materialclass = "dielectric"

    def epsilon(self, xi):
        return 1.0 + 1.5e30 / (1.0e30 + xi**2)
```

The argument `xi` is the imaginary angular frequency in `rad/s`.

## Supported `materialclass` Values

### `"dielectric"`

Use this for insulating or dielectric materials with a finite static permittivity.

Expected behavior:

- `epsilon(0)` is finite
- the TE zero-frequency contribution vanishes in the implementation

### `"drude"`

Use this for metals modeled with dissipation.

In practice, built-in Drude-like materials also provide:

- `wp`: plasma frequency
- `gamma`: damping rate

### `"plasma"`

Use this for dissipationless plasma models.

In practice, these materials provide:

- `wp`: plasma frequency

### `"pec"`

Perfect electric conductor handling is built into the reflection-coefficient logic.

## Notes On Implementation

The package evaluates dielectric response on the imaginary axis, not at real frequencies.

That means user-supplied models should implement `epsilon(xi)` directly for imaginary angular frequency input.

## Layered Materials

To model coatings, pass a list of materials to `matL` or `matR`.

Example:

```python
from califorcia import system
from califorcia.materials import gold, teflon, vacuum

s = system(
    300.0,
    1e-6,
    [teflon, gold],
    gold,
    vacuum,
    deltaL=[50e-9],
)
```

Ordering convention:

- first entry: layer facing the medium
- last entry: substrate half-space

The current implementation supports up to three coating layers on each side.
