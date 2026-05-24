#!/usr/bin/env bash
#
# Set up the Python venv that hosts python-flint (Arb oracle) and
# mpmath (ADR-0035 Tier 2 cross-check oracle) for the pfloat Phase 1
# oracle harness. The venv lives outside the pfloat repo so the LGPL
# FLINT/Arb libraries never enter the shipped Rust crate's link graph
# (ADR-0034 + ADR-0035). mpmath is BSD-licensed; subprocess isolation
# is shared with python-flint for uniformity.
#
# The venv path is `${PFLOAT_ARB_ORACLE_VENV}` if set, otherwise
# `${HOME}/.cache/pfloat-arb-oracle/venv`. Both the Rust ArbOracle
# and MpmathOracle resolve the venv from the same env var / default.
#
# Idempotent: if the venv already has both python-flint and mpmath
# installed the script reports OK and exits. Otherwise it creates
# (or updates) the venv and pip-installs the missing packages.
#
# python-flint is not packaged in nixpkgs (only the C library
# `flint` is, and not the Python binding). mpmath is pure Python and
# trivially pip-installable. The macOS arm64 PyPI wheel for
# python-flint bundles its own flint and arb so no system flint is
# required; on Linux x86_64 the wheel similarly bundles its deps.
# The setup runs without nix-shell.

set -euo pipefail

VENV_PATH="${PFLOAT_ARB_ORACLE_VENV:-${HOME}/.cache/pfloat-arb-oracle/venv}"
PYTHON_BIN="${PFLOAT_ARB_ORACLE_PYTHON:-python3}"

needs_setup=0
if [ ! -d "${VENV_PATH}" ] || [ ! -x "${VENV_PATH}/bin/python3" ]; then
    needs_setup=1
elif ! "${VENV_PATH}/bin/python3" -c "from flint import arb" >/dev/null 2>&1; then
    needs_setup=1
elif ! "${VENV_PATH}/bin/python3" -c "import mpmath" >/dev/null 2>&1; then
    needs_setup=1
fi

if [ "$needs_setup" -eq 0 ]; then
    echo "OK: oracle venv already set up at ${VENV_PATH} (python-flint + mpmath)"
    exit 0
fi

echo "Creating / updating oracle venv at ${VENV_PATH} ..."
mkdir -p "$(dirname "${VENV_PATH}")"
if [ ! -d "${VENV_PATH}" ]; then
    "${PYTHON_BIN}" -m venv "${VENV_PATH}"
fi
"${VENV_PATH}/bin/pip" install --quiet --upgrade pip
"${VENV_PATH}/bin/pip" install --quiet python-flint mpmath

echo "Verifying python-flint install:"
"${VENV_PATH}/bin/python3" -c "
from flint import arb, ctx
ctx.prec = 64
ai, ai_p, bi, bi_p = arb('1.5').airy()
print(f'OK: arb(1.5).airy() = (Ai={ai!s}, Ai\\'={ai_p!s}, Bi={bi!s}, Bi\\'={bi_p!s})')
"

echo "Verifying mpmath install:"
"${VENV_PATH}/bin/python3" -c "
import mpmath
mpmath.mp.prec = 200
print(f'OK: mpmath.airyai(1.5) = {mpmath.airyai(1.5)}')
"
