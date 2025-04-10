# h5-to-pdb
Script for converting ReaDDy output (.h5) into PDB output (.pdb)

# How to use
1. Clone the repository
2. Install the required packages
```bash
python3 -m pip install -editable . # install
python3 -m pip uninstall -u # uninstall
```
3. (optional) make file `mol_settings.json` and write molecule settings
3. Use Python as interpreter (or write in script)
4. Run the script with the following command:
```python3
import h5toPDB as h5pdb
result = h5pdb.h5toPDB('input.h5')
result.generate_pdb() # making pdb file
result.generate_tcl() # making tcl file
result.generate_lammpstrj() # making lammps trajectory file
```
