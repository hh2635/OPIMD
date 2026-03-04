#!/bin/bash

PARAMS="../inputs/parameters.json"

mapfile -t VARS < <(python3 - "$PARAMS" <<'PY'
import json, sys
p = sys.argv[1]
with open(p) as f:
    cfg = json.load(f)

sp = cfg.get("simulation_parameters", {})
pc = cfg.get("physical_constants", {})

pot = sp.get("potential_type", cfg.get("potential_type"))
if pot is None:
    raise SystemExit("Error: potential_type not found (neither in simulation_parameters nor root).")

# beta: root.beta > sp.beta > 1/(kB*T)
beta = cfg.get("beta", sp.get("beta"))
if beta is None:
    T  = sp.get("T")
    kB = pc.get("kB", 1.0)
    if T is None:
        raise SystemExit("Error: neither beta nor T found in config.")
    beta = 1.0 / (kB * T)

print(pot)
print(f"{beta:.2f}")
PY
)

POT="${VARS[0]}"
BETA2="${VARS[1]}"

OUTDIR="../../results/${POT}_beta=${BETA2}"
if [ -d "$OUTDIR" ]; then
  echo "Removing existing directory: $OUTDIR"
  rm -r -- "$OUTDIR"
fi
echo "Creating directory: $OUTDIR"
mkdir -p -- "$OUTDIR"

# sbatch --output="${OUTDIR}/slurm_%j.out" submit_job.slurm
sbatch --job-name="${POT}_beta=${BETA2}" \
       --output="${OUTDIR}/slurm_%j.out" \
       --error="${OUTDIR}/slurm_%j.err" \
       submit_job.slurm

echo "Output is writing to: ${OUTDIR}/slurm_%j.out"