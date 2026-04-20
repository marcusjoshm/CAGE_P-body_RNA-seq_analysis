#!/bin/bash
# Run all Python notebooks sequentially using the cage-pbody kernel
set -e
cd "$(dirname "$0")"

for nb in notebooks/0[1-9]_*.ipynb; do
    echo "=== Running $nb ==="
    jupyter nbconvert --to notebook --execute "$nb" \
        --output "$(basename "$nb")" \
        --ExecutePreprocessor.kernel_name=cage-pbody \
        --ExecutePreprocessor.timeout=600
    echo "=== Completed $nb ==="
done

echo ""
echo "Pipeline complete. Run notebooks/10_plots.Rmd in RStudio for plot reproduction."
