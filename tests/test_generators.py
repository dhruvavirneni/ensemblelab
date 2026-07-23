from ensemblelab import Ensemble, generate


def test_generate_builds_aligned_rdkit_and_ase_conformers() -> None:
    ensemble = generate("CCO", n_confs=3)

    assert isinstance(ensemble, Ensemble)
    assert ensemble.smiles == "CCO"
    assert ensemble.molecule.GetNumConformers() == 3
    assert len(ensemble.conformers) == 3
    assert [conformer.energy for conformer in ensemble.conformers] == [None, None, None]
    assert ensemble.metadata["generator"] == "rdkit.ETKDGv3"
    assert ensemble.metadata["optimization_status"] == "unoptimized"
    assert ensemble.metadata["energy_status"] == "uncomputed"

    for conformer in ensemble.conformers:
        assert conformer.atoms.get_chemical_symbols() == [
            atom.GetSymbol() for atom in ensemble.molecule.GetAtoms()
        ]
        assert len(conformer.atoms) == ensemble.molecule.GetNumAtoms()
        assert conformer.energy is None
        assert ensemble.rdkit_conformer(conformer.id).GetId() == conformer.id


def test_ensemble_from_smiles_delegates_to_generator() -> None:
    ensemble = Ensemble.from_smiles("CCO", n_confs=2)
    assert len(ensemble.conformers) == 2


def test_generate_rejects_invalid_inputs() -> None:
    for invalid_smiles in ("", "not a smiles"):
        try:
            generate(invalid_smiles)
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid SMILES should raise ValueError")

    for n_confs in (0, -1, True, 1.5):
        try:
            generate("CCO", n_confs=n_confs)  # type: ignore[arg-type]
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid conformer count should raise ValueError")
