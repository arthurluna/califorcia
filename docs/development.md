# Development Notes

This page is for contributors and maintainers working on the package itself.

## Local Setup

Create and activate a virtual environment, then install dependencies:

```bash
python -m venv .venv
```

- Windows:

```bash
.venv\Scripts\activate
```

- macOS/Linux:

```bash
source .venv/bin/activate
```

Install dependencies and the package in editable mode:

```bash
pip install -r requirements.txt
pip install -e .
```

## Running Tests

```bash
pytest
```

The current regression tests live in [`tests/test_results.py`](../tests/test_results.py).

They verify:

- the perfect-conductor zero-temperature limit against Casimir's exact result
- zero-frequency high-temperature limits for selected systems

## Source Layout

- [`califorcia/compute.py`](../califorcia/compute.py): main `system` class and observable dispatch
- [`califorcia/plane.py`](../califorcia/plane.py): Fresnel and multilayer reflection coefficients
- [`califorcia/frequency_summation.py`](../califorcia/frequency_summation.py): finite-temperature summation methods
- [`califorcia/interaction.py`](../califorcia/interaction.py): integrands and observable-specific kernels
- [`califorcia/materials/`](../califorcia/materials): built-in material models
- [`examples/`](../examples): usage examples
- [`tests/`](../tests): regression tests

## Documentation Maintenance

When updating the implementation, it is worth keeping these areas synchronized:

- public constructor signature of `system`
- supported `observable` and `fs` values
- list of built-in materials
- coating-layer ordering conventions
- numerical assumptions and units
