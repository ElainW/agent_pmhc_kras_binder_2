#!/bin/bash
# One-time environment fixes for dl_binder_design on this Vast.ai instance.
# See pmhc_design/README.md "Environment notes" for why each step is needed.
set -euo pipefail

DL_BINDER_DESIGN=/workspace/dl_binder_design
PMHC_DESIGN="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# mpnn_fr/dl_interface_design.py expects a ProteinMPNN checkout inside mpnn_fr/.
# Reuse the project's own copy instead of re-cloning.
if [ ! -e "$DL_BINDER_DESIGN/mpnn_fr/ProteinMPNN" ]; then
    ln -s /workspace/ProteinMPNN "$DL_BINDER_DESIGN/mpnn_fr/ProteinMPNN"
fi

# dl_interface_design.py needs torch (ProteinMPNN) AND pyrosetta (FastRelax) in the same
# process; the dl_binder_design venv shipped with pyrosetta+jax+tensorflow but no torch.
/venv/dl_binder_design/bin/pip show torch >/dev/null 2>&1 || \
    /venv/dl_binder_design/bin/pip install torch --index-url https://download.pytorch.org/whl/cu128

# AF2 initial-guess model weights (not present in the base image).
PARAMS_DIR="$DL_BINDER_DESIGN/af2_initial_guess/model_weights/params"
if [ ! -f "$PARAMS_DIR/params_model_1_ptm.npz" ]; then
    mkdir -p "$PARAMS_DIR"
    TMP_TAR="$(mktemp -d)/alphafold_params_2022-12-06.tar"
    wget -q https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar -O "$TMP_TAR"
    tar --extract --file="$TMP_TAR" -C "$PARAMS_DIR" params_model_1_ptm.npz
    rm -f "$TMP_TAR"
fi

# RTX 5090 (Blackwell, sm_120) + vendored 2022-era AlphaFold/jax compatibility.
/venv/af2_binder_design/bin/pip show jax | grep -q "Version: 0.10" || \
    /venv/af2_binder_design/bin/pip install -U "jax[cuda12]"
if ! git -C "$DL_BINDER_DESIGN" apply --reverse --check \
        "$PMHC_DESIGN/patches/dl_binder_design_blackwell_jax_compat.patch" >/dev/null 2>&1; then
    git -C "$DL_BINDER_DESIGN" apply "$PMHC_DESIGN/patches/dl_binder_design_blackwell_jax_compat.patch"
fi

echo "dl_binder_design environment ready."
