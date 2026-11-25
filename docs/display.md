# Display

## Overview
`Ensemble.show()` and `Conformer.show()` print human readable summaries of stored ensemble data. They never optimize conformers, compute descriptors, or otherwise perform analysis. The default ensemble view contains the molecular summary and a compact conformer table; optional sections expose workflow provenance and raw metadata for inspection and debugging.

## Functions
- **Ensemble.show():** prints a default ensemble summary and conformer table.
- **Conformer.show():** prints stored energy, optimization, convergence, and atom-count data for one conformer.

### `Ensemble.show()`
**Overview**
The default display reports SMILES, explicit atom count, net charge (e.g. `Charge: -1`), spin multiplicity, conformer count, whether energy values are available, the stored optimization method, and one row per conformer. When every conformer has an energy in one unit, the table shows relative energy (Delta E) from the lowest energy conformer; this keeps the comparison compact while the absolute energies remain on each `Conformer`.

Set `history=True` to add a workflow history section. The display renders the chronological provenance log in `metadata["history"]`; generation, optimization, and filtering records may contain different fields, and only populated fields are shown. Set `metadata=True` to add a raw key/value metadata table. Sections can be combined, and `conformers=False` suppresses the default conformer table.

**Usage**
```python
from ensemblelab import Ensemble

ensemble = Ensemble.from_smiles("CCO", n_confs=25)
ensemble.show()

ensemble.show(history=True, metadata=True, conformers=False)
```

### `Conformer.show`
**Overview**
Returns the conformer ID, stored energy and unit, optimization method, convergence state, and number of atoms. It does not include a coordinate table; molecular geometry is retained in `conformer.atoms` for programmatic use and future visualization tools.

**Usage**

```python
conformer = ensemble.conformers[0]
conformer.show()
```

---

## Visualization (v0 — untested)

### Overview
Interactive 3D visualization of conformer structures using ASE's visualize infrastructure. This module is in early development (v0) and has not been thoroughly tested across environments and viewer backends.

### `Conformer.view()`
**Status**: v0 — experimental, untested

**Overview**
Display the 3D structure of a single conformer using an interactive molecular viewer.

Supported viewers:
- `"ase"`: ASE's built-in GUI (default; requires X11 or equivalent display server)
- `"ngl"`: NGLView Jupyter widget (requires nglview; suitable for notebook environments)

The viewer is launched as a separate process or widget and does not block notebook execution when using NGLView.

**Parameters**
- `viewer` (str, default `"ase"`): The visualization backend to use.

**Usage**

```python
from ensemblelab import Ensemble

ensemble = Ensemble.from_smiles("CCO", n_confs=5)
conformer = ensemble.conformers[0]

# Display with ASE GUI (interactive, desktop)
conformer.view(viewer="ase")

# Display with NGLView (notebook-friendly)
conformer.view(viewer="ngl")
```

**Known limitations**
- ASE GUI viewer requires a display server and may not work in remote or headless environments.
- NGLView requires explicit installation (`pip install nglview`) and works best in Jupyter notebooks.
- Display server configuration (`DISPLAY` environment variable) may need adjustment on Linux systems.
- Conformer geometry must be populated in `conformer.atoms` (ASE `Atoms` object); empty or invalid structures will raise an error.

**Future development**
- Support for additional viewers (PyMOL, py3Dmol, etc.)
- Ensemble-level visualization with conformer comparison
- Trajectory playback for conformer ensembles
- Custom styling and annotation options
