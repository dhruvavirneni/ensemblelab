"""Generate an unoptimized caffeine conformer ensemble."""

from ensemblelab import generate

CAFFEINE = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"
ensemble = generate(CAFFEINE, n_confs=50)
print(ensemble.metadata)
