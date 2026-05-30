#!/usr/bin/env python3
"""One-command smoke test for the paper's verification scripts."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def run(label: str, args: list[str]) -> None:
    print(f"=== {label} ===")
    subprocess.run([sys.executable, *args], cwd=ROOT, check=True)
    print()


def main() -> int:
    run("theorem certificate", ["certify_theorem.py"])
    run("seed certificates", ["verify_q0.py"])
    run("M regression", ["test_m_value.py"])
    print("All smoke checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
