"""Generate an unoptimized ethanol conformer ensemble."""

from ensemblelab import generate

ensemble = generate("CCO", n_confs=20)
print(ensemble.metadata)
