#!/bin/sh
# Single confirmatory run of the registered analysis (prereg/PREREG.md, Section 8). Not part of the frozen code.
{
echo "registration: doi:10.5281/zenodo.24223802 (created 2026-11-14T15:59:33Z)"
echo "run start (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
python3 -c "import sys, numpy, scipy, pandas, pyarrow; print('python', sys.version.split()[0], 'numpy', numpy.__version__, 'scipy', scipy.__version__, 'pandas', pandas.__version__, 'pyarrow', pyarrow.__version__)"
echo "threads: OMP_NUM_THREADS=$OMP_NUM_THREADS OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS; cpus $(nproc)"
sha256sum screen_lib.py transfer_screen.py crosswalk_types.json banc_888_meta.feather banc_888_edgelist_simple_v3.feather
echo "command: python3 transfer_screen.py banc_888_id banc_888_meta.feather banc_888_edgelist_simple_v3.feather results_banc.json 1 5"
echo "----"
python3 transfer_screen.py banc_888_id banc_888_meta.feather banc_888_edgelist_simple_v3.feather results_banc.json 1 5
echo "----"
echo "exit code: $?"
echo "run end (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > run.log 2>&1
