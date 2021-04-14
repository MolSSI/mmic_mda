[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/MolSSI/mmic_mda/workflows/CI/badge.svg)](https://github.com/MolSSI/mmic_mda/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/MolSSI/mmic_mda/branch/main/graph/badge.svg)](https://codecov.io/gh/MolSSI/mmic_mda/branch/main)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/MolSSI/mmic_mda.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/MolSSI/mmic_mda/context:python)

MDAnalysis translator for MMSchema
==============================
This is part of the [MolSSI](http://molssi.org) Molecular Mechanics Interoperable Components ([MMIC](https://github.com/MolSSI/mmic)) project. This package provides translators between MMSchema and [MDAnalysis](https://github.com/MDAnalysis/mdanalysis).

![image](mmic_mda/data/imgs/component.png)

# API
**mmic_mda** provides 3 classes of translators for: molecules, trajectories, and forcefields.

## Models
```python
# Import models for MDAnalysis
from mmic_mda.models import MdaMol, MdaFF, MdaTraj
 
# Create mda objects from MMSchema models
mda_mol = MdaMol.from_schema(mm_mol)
mda_ff = MdaFF.from_schema(mm_ff)
mda_traj = MdaTraj.from_schema(mm_traj)
 
# Convert mda objects to MMSchema models
mm_mol = mda_mol.to_schema()
mm_ff = mda_ff.to_schema()
mm_traj = mda_traj.to_schema()
```

## Components
Converting from `MDAnalysis.Universe` to `MMElemental.models.Molecule` is achieved with `MdaToMolComponent`:

```python
 
# Create input for converting Mda to MMSchema molecule
trans_in = {
        "data_object": MDAnalysis.Universe,
        "data_units": mmic_mda.units
}
 
# Run conversion
trans_out = MdaToMolComponent.compute(trans_in)

mm_mol = trans_out.schema_object -> MMElemental.models.Molecule
```

Converting from `MMElemental.models.Molecule` to `MDAnalysis.Universe` is achieved with `MolToMdaComponent`:
 
# Create input for converting MMSchema to Mda molecule
trans_in = {
        "schema_object": MDAnalysis.Universe,
        "engine_version": '>=1.0.0' # reject conversion for MDAv < 1.0.0
)
 
# Run conversion
trans_out = MolToMdaComponent.compute(trans_in)

mda_uni = trans_out.data_object -> MDAnalysis.Universe
```

### Copyright
Copyright (c) 2021, MolSSI


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
