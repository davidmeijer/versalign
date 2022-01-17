#!/usr/bin/env python3
"""
Author: David Meijer
Description: test script.
"""
import sys
from pathlib import Path
path_root = Path(__file__).parents[1]
sys.path.append(str(path_root))

from monomer_aligner import run_aligner


if __name__ == '__main__':
    run_aligner()
