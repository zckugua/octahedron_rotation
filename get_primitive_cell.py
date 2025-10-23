from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

structure = Poscar.from_file("POSCAR_AVG").structure
analyzer = SpacegroupAnalyzer(structure, symprec=0.03)
primitive_structure = analyzer.get_primitive_standard_structure()
conventional_structure = analyzer.get_conventional_standard_structure()
primitive_structure.to(fmt="poscar", filename="POSCAR_primitive")
conventional_structure.to(fmt="poscar", filename="POSCAR_conventional")

