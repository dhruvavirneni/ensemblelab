import ensemblelab.generators as generators


def test_generate_import_works_without_eager_optimizer_imports():
    ensemble = generators.generate("CCO", n_confs=2)
    assert len(ensemble.conformers) == 2
    assert ensemble.smiles == "CCO"
