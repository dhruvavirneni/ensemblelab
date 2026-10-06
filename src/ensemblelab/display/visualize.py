from ..generators import Ensemble, Conformer
import ase.visualize

# include link to ASE view() in documentation...outline viewing methods in docs though
# viewer options include ase.gui and ngl for now
# figure out other parameters for ase view
def view(conformer: Conformer, viewer: str = "ase") -> None:
    """3D visualization of inputted conformer; built from ASE's visualize library"""
    VIEWER_LIST = ["ase", "ngl"]
    # error handling
    if len(conformer.atoms) == 0:
        raise ValueError("Conformer does not contain calculated ASE Atoms")

    if viewer not in VIEWER_LIST:
        raise ValueError("Invalid viewer type.")

    atoms = conformer.atoms
    ase.visualize.view(atoms=atoms, viewer=viewer)