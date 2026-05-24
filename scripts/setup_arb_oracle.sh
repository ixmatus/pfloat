#!/usr/bin/env bash
#
# Set up the Python venv that hosts python-flint for the pfloat
# Arb oracle backend. The venv lives outside the pfloat repo so the
# LGPL FLINT/Arb libraries never enter the shipped Rust crate's
# link graph (ADR-0034).
#
# The venv path is `${PFLOAT_ARB_ORACLE_VENV}` if set, otherwise
# `${HOME}/.cache/pfloat-arb-oracle/venv`. The Rust ArbOracle
# resolves the venv from the same env var / default.
#
# Idempotent: if the venv already has python-flint installed the
# script just reports OK and exits. Otherwise it creates the venv
# (via the system `python3`) and pip-installs python-flint.
#
# python-flint is not packaged in nixpkgs (only the C library
# `flint` is, and not the Python binding). The macOS arm64 PyPI
# wheel bundles its own flint and arb so no system flint is
# required; on Linux x86_64 the wheel similarly bundles its deps.
# The setup runs without nix-shell.

set -euo pipefail

VENV_PATH="${PFLOAT_ARB_ORACLE_VENV:-${HOME}/.cache/pfloat-arb-oracle/venv}"
PYTHON_BIN="${PFLOAT_ARB_ORACLE_PYTHON:-python3}"

if [ -d "${VENV_PATH}" ] \
   && [ -x "${VENV_PATH}/bin/python3" ] \
   && "${VENV_PATH}/bin/python3" -c "from flint import arb" >/dev/null 2>&1; then
    echo "OK: Arb oracle venv already set up at ${VENV_PATH}"
    exit 0
fi

echo "Creating Arb oracle venv at ${VENV_PATH} ..."
mkdir -p "$(dirname "${VENV_PATH}")"
"${PYTHON_BIN}" -m venv "${VENV_PATH}"
"${VENV_PATH}/bin/pip" install --quiet --upgrade pip
"${VENV_PATH}/bin/pip" install --quiet python-flint

echo "Verifying python-flint install:"
"${VENV_PATH}/bin/python3" -c "
from flint import arb, ctx
ctx.prec = 64
ai, ai_p, bi, bi_p = arb('1.5').airy()
print(f'OK: arb(1.5).airy() = (Ai={ai!s}, Ai\\'={ai_p!s}, Bi={bi!s}, Bi\\'={bi_p!s})')
"
