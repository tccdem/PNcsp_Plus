#!/usr/bin/env python3
import os
import sys
import json
import glob
import logging
import argparse
import subprocess
import warnings
from typing import Dict, Any, Tuple, List

import numpy as np
from ase.io import read, write
from ase.data import atomic_masses, atomic_numbers

warnings.filterwarnings("ignore")
from pymatgen.core import Structure, Element

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

OTFG_STRINGS = {
    'Ag': '3|1.6|16.5|18|19|40U:50UU:41UU:42UU(qc=7)[]',
    'Al': '3|2|3.5|5|7|30UU:31UU:32UU[]',
    'As': '3|2|11|13|15|40UU:41UU:32UU(qc=6)[]',
    'Au': '3|2.4|9.2|11|11.8|50U:60UU:51UU:52UU:43UU(qc=6)[]',
    'B': '2|1.2|7.4|10|12.9|20UU:21UU(qc=8)[]',
    'Ba': '2|2|4|5.5|7|50U:60UU:51UU(qc=5.5)[]',
    'Be': '2|1|13.5|17.6|21|10U:20UU:21UU(qc=7)[]',
    'Bi': '3|2.3|8|9.5|11|60UU:61UU:52UU[]',
    'Br': '2|2|4|5.5|7|40UU:41UU(qc=6)[]',
    'C': '2|1.4|11|12|14.7|20UU:21UU(qc=7)[]',
    'Ca': '3|2|8|10|12|30U:40U:31UU:32UU[]',
    'Cd': '1|2.2|8.7|9.5|11|50U+0U+0.1:42UU(qc=5,q0=4)[]',
    'Ce': '2|2.1|9|12|14|50U:60UU:51UU:52L:43UU(qc=6)[]',
    'Cl': '3|1.8|6|7|8.5|30UU:31UU:32UU[]',
    'Co': '3|2.2|2|1|9|12|14|40UU:32UU:41UU(qc=5.5)[]',
    'Cr': '3|2|2|0.7|11|12|14|30U:40UU:31UU:32UU(qc=6)[]',
    'Cs': '2|2.2|3|4|5|50U:60UU:51UU(qc=5)[]',
    'Cu': '3|2.2|2|1|10|13|15|40UU:41UU:32UU(qc=6)[]',
    'Dy': '2|1.9|14|16|18|50U:60UU:51UU:43UU(qc=6.5)[]',
    'Er': '2|2.1|11|13|14|50U:60UU:51UU:43UU{6s0.5}(qc=6)[]',
    'Eu': '2|1.9|12|14|16|50U:60UU:51UU:43UU(qc=6)[]',
    'F': '2|1.2|10.5|13|14.5|20UU:21UU(qc=7)[]',
    'Fe': '3|2.2|2|1|9|11|13|40UU:41UU:32UU(qc=5.5)[]',
    'Ga': '3|2|13|14|15|40UU:41UU:32UU(qc=6)[]',
    'Gd': '2|1.8|13|14.7|16.5|50U:60UU:51UU:43UU{5d0,4f8}(qc=6.5)[]',
    'Ge': '3|2|11|13|14.5|40UU:41UU:32UU(qc=6)[]',
    'H': '1|0.6|9|10|15|10UU(qc=8)[]',
    'Hf': '3|2.1|12|14|16|50U:52UU:60UU:51UU:43UU(qc=6)[]',
    'Hg': '3|2.2|8.5|10|11.4|60UU:61UU:50U:51U:52UU[]',
    'Ho': '2|1.9|14|16|18|50U:60UU:51UU:43UU(qc=6.5)[]',
    'I': '2|2|4|6|8|50UU:51UU[]',
    'In': '3|2.3|9|10.5|12|50UU:51UU:42UU[]',
    'Ir': '3|2.4|8.8|10|12|50U:60UU:51UU:52UU:43UU(qc=6)[]',
    'K': '2|1.5|10|13|16|30U:40UU:31UU(qc=6)[]',
    'La': '2|2.3|5|6|8.8|50U:60UU:51UU:52UU(qc=4.5)[]',
    'Li': '1|1|12|15|20|10U:20UU(qc=7)[]',
    'Lu': '2|2.1|11|14|15|50U:60UU:51UU:52UU:43UU(qc=6)[]',
    'Mg': '3|1.8|10|14|18|20U:30UU:21UU:32UU[]',
    'Mn': '3|2.2|2|0.7|9|10.7|12.1|40UU:32UU:41UU(qc=5.5)[]',
    'Mo': '3|1.6|11|14|15.5|40U:50UU:41UU:42UU(qc=6)[]',
    'N': '2|1.1|12|16|19|20UU:21UU(qc=7)[]',
    'Na': '2|1.3|15|19|22|20U:30UU:21UU(qc=7)[]',
    'Nb': '3|1.6|10|12.9|14.7|40U:50UU:41UU:42UU(qc=6)[]',
    'Nd': '2|2|13|15|16|50U:60UU:51UU:43UU(qc=6)[]',
    'Ni': '3|2.2|2|1|10|12.1|14.7|40UU:41UU:32UU(qc=6)[]',
    'O': '2|1.1|14|18|21|20UU:21UU(qc=7)[]',
    'Os': '3|2.4|10|12.1|13.6|50U:60UU:51UU:52UU:43UU(qc=6)[]',
    'P': '3|1.8|2|4|6|30UU:31UU:32UU[]',
    'Pb': '3|2.2|2.5|1.8|8.1|9|10.5|60UU:61UU:52UU:51U:50U[]',
    'Pd': '3|1.6|15|18.8|20.6|40U:50UU:41UU:42UU{5s0.5}(qc=7)[]',
    'Pm': '2|2|11|14|15|50U:60UU:51UU:43UU(qc=6)[]',
    'Po': '3|2.3|8.1|10|11.5|60UU:61UU:52UU[]',
    'Pr': '2|2.1|9|11.8|14.7|50U:60UU:51UU:43UU(qc=6)[]',
    'Pt': '3|2.4|8.5|10|11.8|50U:60UU:51UU:52UU:43UU(qc=6)[]',
    'Rb': '2|2.1|6|7.4|9.2|40U:50UU:41UU[]',
    'Re': '3|2.4|8.8|11|12.8|50U:60UU:51UU:52UU:43UU(qc=6)[]',
    'Rh': '3|1.6|14.3|18.4|20.2|50U:40U:41UU:42UU(qc=7)[]',
    'Ru': '3|1.6|12|14.7|16.5|40U:50UU:41UU:42UU(qc=6)[]',
    'S': '3|1.8|5|6.5|8|30UU:31UU:32UU[]',
    'Sb': '3|2.2|8|9|11|50UU:51UU:42UU[]',
    'Sc': '3|1.8|11|13.5|15.1|30U:40UU:31UU:32UU(qc=6.5)[]',
    'Se': '3|2|10|13|14.5|40UU:41UU:32UU(qc=6)[]',
    'Si': '3|1.8|3.6|5|7|30UU:31UU:32UU[]',
    'Sm': '2|1.9|12|15|17|50U:60UU:51UU:43UU(qc=6)[]',
    'Sn': '3|2.2|8.5|10|11.5|50UU:51UU:42UU[]',
    'Sr': '3|2|7|8.5|10|40U:50UU:41UU:42UU[]',
    'Ta': '3|2.4|9|11.4|13.2|50U:52UU:60UU:51UU:43UU(qc=6)[]',
    'Tb': '2|1.8|16|18|20.2|50U:60UU:51UU:43UU(qc=6.5)[]',
    'Te': '2|2.2|3.5|5.5|7.5|50UU:51UU[]',
    'Ti': '3|1.8|10|12|14|30U:40UU:31UU:32UU(qc=5.5)[]',
    'Tl': '3|2.2|10|11.5|12.5|60UU:61UU:50U:51U:52UU[]',
    'Tm': '2|2.1|12|14|17|50U:60UU:51UU:43UU{4f12}(qc=6)[]',
    'V': '3|2|2|1|10|12|14|30U:40UU:31UU:32UU(qc=6)[]',
    'W': '3|2.4|9.5|11|13|50U:60UU:51UU:52UU:43UU(qc=6)[]',
    'Y': '3|2|6.5|8.5|10|40U:50UU:41UU:42UU[]',
    'Yb': '2|2.1|11.4|13.6|16|50U:60UU:51UU:43UU{4f13}(qc=6)[]',
    'Zn': '3|2|2|1|11|12.8|14.5|40UU:41UU:32UU(qc=6)[]',
    'Zr': '3|2.1|7|8.5|10|40U:50UU:41UU:42UU[]'
}

INTERNAL_CUTOFFS = {
    'Ac': {'coarse': 100.0, 'medium': 200.0, 'fine': 300.0},
    'Ag': {'coarse': 300.0, 'medium': 400.0, 'fine': 500.0},
    'Al': {'coarse': 100.0, 'medium': 150.0, 'fine': 200.0},
    'Ar': {'coarse': 200.0, 'medium': 275.0, 'fine': 350.0},
    'As': {'coarse': 200.0, 'medium': 250.0, 'fine': 300.0},
    'At': {'coarse': 370.0, 'medium': 430.0, 'fine': 500.0},
    'Au': {'coarse': 300.0, 'medium': 400.0, 'fine': 500.0},
    'B': {'coarse': 250.0, 'medium': 320.0, 'fine': 400.0},
    'Ba': {'coarse': 200.0, 'medium': 300.0, 'fine': 400.0},
    'Be': {'coarse': 200.0, 'medium': 300.0, 'fine': 400.0},
    'Bi': {'coarse': 340.0, 'medium': 400.0, 'fine': 500.0},
    'Br': {'coarse': 120.0, 'medium': 175.0, 'fine': 250.0},
    'C': {'coarse': 300.0, 'medium': 470.0, 'fine': 600.0},
    'Ca': {'coarse': 300.0, 'medium': 400.0, 'fine': 500.0},
    'Cd': {'coarse': 330.0, 'medium': 450.0, 'fine': 600.0},
    'Ce': {'coarse': 400.0, 'medium': 550.0, 'fine': 700.0},
    'Cl': {'coarse': 160.0, 'medium': 300.0, 'fine': 450.0},
    'Co': {'coarse': 400.0, 'medium': 450.0, 'fine': 550.0},
    'Cr': {'coarse': 400.0, 'medium': 500.0, 'fine': 700.0},
    'Cs': {'coarse': 200.0, 'medium': 300.0, 'fine': 400.0},
    'Cu': {'coarse': 320.0, 'medium': 450.0, 'fine': 600.0},
    'F': {'coarse': 630.0, 'medium': 780.0, 'fine': 850.0},
    'Fe': {'coarse': 340.0, 'medium': 400.0, 'fine': 500.0},
    'Fr': {'coarse': 200.0, 'medium': 300.0, 'fine': 400.0},
    'Ga': {'coarse': 520.0, 'medium': 590.0, 'fine': 800.0},
    'Gd': {'coarse': 700.0, 'medium': 740.0, 'fine': 820.0},
    'Ge': {'coarse': 120.0, 'medium': 180.0, 'fine': 300.0},
    'H': {'coarse': 250.0, 'medium': 300.0, 'fine': 350.0},
    'He': {'coarse': 620.0, 'medium': 1100.0, 'fine': 1400.0},
    'Hf': {'coarse': 200.0, 'medium': 250.0, 'fine': 300.0},
    'Hg': {'coarse': 340.0, 'medium': 400.0, 'fine': 450.0},
    'I': {'coarse': 200.0, 'medium': 300.0, 'fine': 400.0},
    'In': {'coarse': 50.0, 'medium': 100.0, 'fine': 200.0},
    'Ir': {'coarse': 340.0, 'medium': 400.0, 'fine': 500.0},
    'K': {'coarse': 270.0, 'medium': 340.0, 'fine': 400.0},
    'Kr': {'coarse': 150.0, 'medium': 200.0, 'fine': 300.0},
    'La': {'coarse': 100.0, 'medium': 200.0, 'fine': 300.0},
    'Li': {'coarse': 150.0, 'medium': 300.0, 'fine': 450.0},
    'Lu': {'coarse': 500.0, 'medium': 750.0, 'fine': 900.0},
    'Mg': {'coarse': 100.0, 'medium': 200.0, 'fine': 300.0},
    'Mn': {'coarse': 300.0, 'medium': 350.0, 'fine': 400.0},
    'Mo': {'coarse': 300.0, 'medium': 350.0, 'fine': 400.0},
    'N': {'coarse': 400.0, 'medium': 550.0, 'fine': 700.0},
    'Na': {'coarse': 400.0, 'medium': 550.0, 'fine': 700.0},
    'Nb': {'coarse': 220.0, 'medium': 280.0, 'fine': 350.0},
    'Ne': {'coarse': 400.0, 'medium': 450.0, 'fine': 500.0},
    'Ni': {'coarse': 500.0, 'medium': 550.0, 'fine': 600.0},
    'O': {'coarse': 330.0, 'medium': 450.0, 'fine': 500.0},
    'Os': {'coarse': 300.0, 'medium': 350.0, 'fine': 400.0},
    'P': {'coarse': 400.0, 'medium': 450.0, 'fine': 500.0},
    'Pb': {'coarse': 400.0, 'medium': 500.0, 'fine': 700.0},
    'Pd': {'coarse': 350.0, 'medium': 400.0, 'fine': 450.0},
    'Po': {'coarse': 360.0, 'medium': 420.0, 'fine': 500.0},
    'Pt': {'coarse': 350.0, 'medium': 400.0, 'fine': 450.0},
    'Ra': {'coarse': 200.0, 'medium': 300.0, 'fine': 400.0},
    'Rb': {'coarse': 150.0, 'medium': 200.0, 'fine': 300.0},
    'Re': {'coarse': 300.0, 'medium': 350.0, 'fine': 400.0},
    'Rh': {'coarse': 320.0, 'medium': 360.0, 'fine': 400.0},
    'Rn': {'coarse': 350.0, 'medium': 400.0, 'fine': 450.0},
    'Ru': {'coarse': 300.0, 'medium': 350.0, 'fine': 400.0},
    'S': {'coarse': 300.0, 'medium': 400.0, 'fine': 450.0},
    'Sb': {'coarse': 150.0, 'medium': 200.0, 'fine': 300.0},
    'Sc': {'coarse': 250.0, 'medium': 300.0, 'fine': 400.0},
    'Se': {'coarse': 250.0, 'medium': 350.0, 'fine': 400.0},
    'Si': {'coarse': 125.0, 'medium': 200.0, 'fine': 300.0},
    'Sm': {'coarse': 400.0, 'medium': 570.0, 'fine': 700.0},
    'Sn': {'coarse': 100.0, 'medium': 150.0, 'fine': 200.0},
    'Sr': {'coarse': 170.0, 'medium': 400.0, 'fine': 600.0},
    'Ta': {'coarse': 250.0, 'medium': 300.0, 'fine': 350.0},
    'Tc': {'coarse': 300.0, 'medium': 350.0, 'fine': 400.0},
    'Te': {'coarse': 100.0, 'medium': 150.0, 'fine': 250.0},
    'Ti': {'coarse': 300.0, 'medium': 350.0, 'fine': 400.0},
    'Tl': {'coarse': 350.0, 'medium': 400.0, 'fine': 450.0},
    'V': {'coarse': 290.0, 'medium': 350.0, 'fine': 400.0},
    'W': {'coarse': 250.0, 'medium': 300.0, 'fine': 450.0},
    'Xe': {'coarse': 250.0, 'medium': 300.0, 'fine': 350.0},
    'Y': {'coarse': 200.0, 'medium': 250.0, 'fine': 300.0},
    'Zn': {'coarse': 380.0, 'medium': 460.0, 'fine': 600.0},
    'Zr': {'coarse': 220.0, 'medium': 280.0, 'fine': 360.0}
}

# ==============================================================================
# ENTERPRISE OTFG CUTOFFS (HARTREES) - Active Cutoff Engine
# ==============================================================================
OTFG_CUTOFFS_HA = {
    'Ac': {'coarse': 5.0, 'medium': 6.0, 'fine': 7.0},
    'Ag': {'coarse': 16.5, 'medium': 18.0, 'fine': 19.0},
    'Al': {'coarse': 3.5, 'medium': 5.0, 'fine': 7.0},
    'Am': {'coarse': 13.0, 'medium': 14.0, 'fine': 17.0},
    'Ar': {'coarse': 10.0, 'medium': 12.0, 'fine': 15.0},
    'As': {'coarse': 11.0, 'medium': 13.0, 'fine': 15.0},
    'At': {'coarse': 9.2, 'medium': 10.5, 'fine': 11.5},
    'Au': {'coarse': 9.2, 'medium': 11.0, 'fine': 11.8},
    'B': {'coarse': 7.4, 'medium': 10.0, 'fine': 12.9},
    'Ba': {'coarse': 4.0, 'medium': 5.5, 'fine': 7.0},
    'Be': {'coarse': 13.5, 'medium': 17.6, 'fine': 21.0},
    'Bi': {'coarse': 8.0, 'medium': 9.5, 'fine': 11.0},
    'Bk': {'coarse': 11.4, 'medium': 13.6, 'fine': 16.0},
    'Br': {'coarse': 4.0, 'medium': 5.5, 'fine': 7.0},
    'C': {'coarse': 11.0, 'medium': 12.0, 'fine': 14.7},
    'Ca': {'coarse': 8.0, 'medium': 10.0, 'fine': 12.0},
    'Cd': {'coarse': 8.7, 'medium': 9.5, 'fine': 11.0},
    'Ce': {'coarse': 9.0, 'medium': 12.0, 'fine': 14.0},
    'Cf': {'coarse': 10.0, 'medium': 11.0, 'fine': 12.0},
    'Cl': {'coarse': 6.0, 'medium': 7.0, 'fine': 8.5},
    'Cm': {'coarse': 12.5, 'medium': 14.5, 'fine': 16.5},
    'Co': {'coarse': 9.0, 'medium': 12.0, 'fine': 14.0},
    'Cr': {'coarse': 11.0, 'medium': 12.0, 'fine': 14.0},
    'Cs': {'coarse': 3.0, 'medium': 4.0, 'fine': 5.0},
    'Cu': {'coarse': 10.0, 'medium': 13.0, 'fine': 15.0},
    'Dy': {'coarse': 14.0, 'medium': 16.0, 'fine': 18.0},
    'Er': {'coarse': 11.0, 'medium': 13.0, 'fine': 14.0},
    'Es': {'coarse': 12.0, 'medium': 14.0, 'fine': 16.0},
    'Eu': {'coarse': 12.0, 'medium': 14.0, 'fine': 16.0},
    'F': {'coarse': 10.5, 'medium': 13.0, 'fine': 14.5},
    'Fe': {'coarse': 9.0, 'medium': 11.0, 'fine': 13.0},
    'Fm': {'coarse': 15.0, 'medium': 17.0, 'fine': 18.0},
    'Fr': {'coarse': 5.0, 'medium': 7.5, 'fine': 9.0},
    'Ga': {'coarse': 13.0, 'medium': 14.0, 'fine': 15.0},
    'Gd': {'coarse': 13.0, 'medium': 14.7, 'fine': 16.5},
    'Ge': {'coarse': 11.0, 'medium': 13.0, 'fine': 14.5},
    'H': {'coarse': 9.0, 'medium': 10.0, 'fine': 15.0},
    'He': {'coarse': 9.0, 'medium': 13.0, 'fine': 18.0},
    'Hf': {'coarse': 12.0, 'medium': 14.0, 'fine': 16.0},
    'Hg': {'coarse': 8.5, 'medium': 10.0, 'fine': 11.4},
    'Ho': {'coarse': 14.0, 'medium': 16.0, 'fine': 18.0},
    'I': {'coarse': 4.0, 'medium': 6.0, 'fine': 8.0},
    'In': {'coarse': 9.0, 'medium': 10.5, 'fine': 12.0},
    'Ir': {'coarse': 8.8, 'medium': 10.0, 'fine': 12.0},
    'K': {'coarse': 10.0, 'medium': 13.0, 'fine': 16.0},
    'Kr': {'coarse': 4.0, 'medium': 6.0, 'fine': 8.0},
    'La': {'coarse': 5.0, 'medium': 6.0, 'fine': 8.8},
    'Li': {'coarse': 12.0, 'medium': 15.0, 'fine': 20.0},
    'Lr': {'coarse': 8.0, 'medium': 10.0, 'fine': 12.0},
    'Lu': {'coarse': 11.0, 'medium': 14.0, 'fine': 15.0},
    'Md': {'coarse': 13.0, 'medium': 15.0, 'fine': 16.0},
    'Mg': {'coarse': 10.0, 'medium': 14.0, 'fine': 18.0},
    'Mn': {'coarse': 9.0, 'medium': 10.7, 'fine': 12.1},
    'Mo': {'coarse': 11.0, 'medium': 14.0, 'fine': 15.5},
    'N': {'coarse': 12.0, 'medium': 16.0, 'fine': 19.0},
    'Na': {'coarse': 15.0, 'medium': 19.0, 'fine': 22.0},
    'Nb': {'coarse': 10.0, 'medium': 12.9, 'fine': 14.7},
    'Nd': {'coarse': 13.0, 'medium': 15.0, 'fine': 16.0},
    'Ne': {'coarse': 10.0, 'medium': 15.0, 'fine': 20.0},
    'Ni': {'coarse': 10.0, 'medium': 12.1, 'fine': 14.7},
    'No': {'coarse': 8.0, 'medium': 10.0, 'fine': 12.0},
    'Np': {'coarse': 13.2, 'medium': 15.4, 'fine': 16.5},
    'O': {'coarse': 14.0, 'medium': 18.0, 'fine': 21.0},
    'Os': {'coarse': 10.0, 'medium': 12.1, 'fine': 13.6},
    'P': {'coarse': 2.0, 'medium': 4.0, 'fine': 6.0},
    'Pa': {'coarse': 10.0, 'medium': 12.5, 'fine': 14.0},
    'Pb': {'coarse': 8.1, 'medium': 9.0, 'fine': 10.5},
    'Pd': {'coarse': 15.0, 'medium': 18.8, 'fine': 20.6},
    'Pm': {'coarse': 11.0, 'medium': 14.0, 'fine': 15.0},
    'Po': {'coarse': 8.1, 'medium': 10.0, 'fine': 11.5},
    'Pr': {'coarse': 9.0, 'medium': 11.8, 'fine': 14.7},
    'Pt': {'coarse': 8.5, 'medium': 10.0, 'fine': 11.8},
    'Pu': {'coarse': 12.0, 'medium': 13.0, 'fine': 14.0},
    'Ra': {'coarse': 4.0, 'medium': 7.0, 'fine': 10.0},
    'Rb': {'coarse': 6.0, 'medium': 7.4, 'fine': 9.2},
    'Re': {'coarse': 8.8, 'medium': 11.0, 'fine': 12.8},
    'Rh': {'coarse': 14.3, 'medium': 18.4, 'fine': 20.2},
    'Rn': {'coarse': 8.1, 'medium': 10.0, 'fine': 11.0},
    'Ru': {'coarse': 12.0, 'medium': 14.7, 'fine': 16.5},
    'S': {'coarse': 5.0, 'medium': 6.5, 'fine': 8.0},
    'Sb': {'coarse': 8.0, 'medium': 9.0, 'fine': 11.0},
    'Sc': {'coarse': 11.0, 'medium': 13.5, 'fine': 15.1},
    'Se': {'coarse': 10.0, 'medium': 13.0, 'fine': 14.5},
    'Si': {'coarse': 3.6, 'medium': 5.0, 'fine': 7.0},
    'Sm': {'coarse': 12.0, 'medium': 15.0, 'fine': 17.0},
    'Sn': {'coarse': 8.5, 'medium': 10.0, 'fine': 11.5},
    'Sr': {'coarse': 7.0, 'medium': 8.5, 'fine': 10.0},
    'Ta': {'coarse': 9.0, 'medium': 11.4, 'fine': 13.2},
    'Tb': {'coarse': 16.0, 'medium': 18.0, 'fine': 20.2},
    'Tc': {'coarse': 11.0, 'medium': 13.6, 'fine': 15.4},
    'Te': {'coarse': 3.5, 'medium': 5.5, 'fine': 7.5},
    'Th': {'coarse': 6.0, 'medium': 12.0, 'fine': 13.0},
    'Ti': {'coarse': 10.0, 'medium': 12.0, 'fine': 14.0},
    'Tl': {'coarse': 10.0, 'medium': 11.5, 'fine': 12.5},
    'Tm': {'coarse': 12.0, 'medium': 14.0, 'fine': 17.0},
    'U': {'coarse': 11.5, 'medium': 14.0, 'fine': 16.0},
    'V': {'coarse': 10.0, 'medium': 12.0, 'fine': 14.0},
    'W': {'coarse': 9.5, 'medium': 11.0, 'fine': 13.0},
    'Xe': {'coarse': 3.0, 'medium': 5.0, 'fine': 6.0},
    'Y': {'coarse': 6.5, 'medium': 8.5, 'fine': 10.0},
    'Yb': {'coarse': 11.4, 'medium': 13.6, 'fine': 16.0},
    'Zn': {'coarse': 11.0, 'medium': 12.8, 'fine': 14.5},
    'Zr': {'coarse': 7.0, 'medium': 8.5, 'fine': 10.0}
}

def load_config(config_path: str) -> Dict:
    with open(config_path, 'r') as f:
        return json.load(f)

def analyze_structure_properties(cif_path: str, config: Dict) -> Dict[str, Any]:
    try:
        struct = Structure.from_file(cif_path)
        comp = struct.composition
        exact_z = int(comp.get_reduced_composition_and_factor()[1])
    except Exception as e:
        logger.error(f"Pymatgen parse failed for {os.path.basename(cif_path)}. Falling back to defaults.")
        return {"is_metal": True, "modulus_est": config["heuristics"].get("default_modulus_gpa", 50.0), "charge": 0.0, "z_factor": 1, "formula": "Unknown", "complexity": 1.0, "packing_fraction": "Unknown"}

    is_metal = any(e.is_transition_metal or e.is_actinoid or e.is_lanthanoid or e.is_alkali or e.is_alkaline for e in comp.elements)

    complexity = 1.0
    if len(struct) > 30:
        complexity = 2.5 

    packing_fraction = 0.5 

    elemental_modulus = sum((amt / comp.num_atoms) * (Element(sym).bulk_modulus or config["heuristics"].get("default_modulus_gpa", 50.0)) for sym, amt in comp.get_el_amt_dict().items())
    scaled_modulus = elemental_modulus * (packing_fraction / 0.74)


    return {
        "is_metal": is_metal,
        "modulus_est": round(max(20.0, min(500.0, scaled_modulus)), 1),
        "charge": struct.charge,
        "z_factor": exact_z,
        "formula": comp.reduced_formula,
        "complexity": round(complexity, 3),
        "packing_fraction": round(packing_fraction, 3)
    }

def get_hunds_rule_unpaired_electrons(sym: str) -> float:
    try:
        el = Element(sym)
        if not (el.is_transition_metal or el.is_lanthanoid or el.is_actinoid): return 0.0
        unpaired = 0.0
        for n, orbital, count in el.full_electronic_structure:
            if orbital == 'd': unpaired += min(count, 10 - count)
            elif orbital == 'f': unpaired += min(count, 14 - count)
        return float(unpaired)
    except Exception: return 0.0

def inject_quantum_spins(cell_text: str) -> Tuple[str, float]:
    total_spin = 0.0
    out, in_block = [], False
    for line in cell_text.splitlines():
        l = line.strip().lower()
        if l.startswith("%block positions_frac"): in_block = True; out.append(line); continue
        if in_block and l.startswith("%endblock positions_frac"): in_block = False; out.append(line); continue
        
        if in_block and line.strip() and not line.strip().startswith(("#", "!")):
            parts = line.split()
            if len(parts) >= 4:
                sym = "".join(c for c in parts[0] if c.isalpha()).capitalize()
                spin_val = get_hunds_rule_unpaired_electrons(sym)
                if spin_val > 0:
                    total_spin += spin_val
                    if "SPIN=" not in line.upper(): line = f"{line} SPIN={spin_val:.5f}"
        out.append(line)
    return "\n".join(out) + "\n", total_spin

def process_one_cif(cif_path: str, config: Dict) -> None:
    dirs = config.get("directories", {})
    seed = os.path.splitext(os.path.basename(cif_path))[0]
    work_dir = os.path.join(os.path.expanduser(dirs.get("output_castep", "./")), seed)
    os.makedirs(work_dir, exist_ok=True)
    
    cell_filename = f"{seed}.cell"
    cell_path = os.path.join(work_dir, cell_filename)
    abs_cif = os.path.abspath(cif_path)
    
    try:
        subprocess.run(["cif2cell", "-f", abs_cif, "-p", "castep", "-o", cell_filename], cwd=work_dir, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        logger.warning(f"  [{seed}] cif2cell native parsing failed (likely space group syntax). Initiating ASE sanitization fallback.")
        clean_cif = os.path.join(work_dir, f"{seed}_sanitized.cif")
        try:
            temp_atoms = read(abs_cif)
            write(clean_cif, temp_atoms, format="cif")
            subprocess.run(["cif2cell", "-f", clean_cif, "-p", "castep", "-o", cell_filename], cwd=work_dir, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if os.path.exists(clean_cif):
                os.remove(clean_cif)
        except Exception as fallback_err:
            logger.error(f"  [{seed}] Fallback sanitization failed: {fallback_err}")
            return

    props = analyze_structure_properties(cif_path, config)
    
    try:
        atoms = read(cif_path)
    except Exception as e:
        logger.error(f"ASE failed to read {seed}. Skipping. Error: {e}")
        return

    symbols = list(dict.fromkeys(atoms.get_chemical_symbols()))
    if props["formula"] == "Unknown": props["formula"] = "".join(symbols)

    logger.info(f"\n--- Processing: {seed} ---")
    
    total_val = sum(Element(a.symbol).group for a in atoms)
    calc_nextra = max(10, int(total_val * 0.1))
    
    kpt_spacing = config.get("heuristics", {}).get("kpoint_spacing", 0.03)
    cutoff_quality = config.get("heuristics", {}).get("cutoff_quality", "ultrafine").lower()
    
    if cutoff_quality not in ["coarse", "medium", "fine", "ultrafine"]:
        logger.warning(f"Invalid cutoff_quality '{cutoff_quality}'. Defaulting to 'ultrafine'.")
        cutoff_quality = "ultrafine"

    max_base_cutoff = 0.0
    HARTREE_TO_EV = 27.211386
    
    dict_quality = 'fine' if cutoff_quality == 'ultrafine' else cutoff_quality
    
    for sym in symbols:
        el_cutoffs = OTFG_CUTOFFS_HA.get(sym, {'coarse': 12.0, 'medium': 15.0, 'fine': 18.0})
        val_ha = el_cutoffs.get(dict_quality, el_cutoffs['fine'])
        val_ev = float(val_ha) * HARTREE_TO_EV
        logger.info(f"  [{sym}] Retrieved Base {dict_quality.upper()} Cutoff: {val_ha} Ha ({val_ev:.1f} eV)")
        max_base_cutoff = max(max_base_cutoff, val_ev)

    if cutoff_quality == 'ultrafine':
        calc_cutoff = np.ceil((max_base_cutoff + 60.0) / 10.0) * 10
        cutoff_log_msg = f"Dynamic Buffered Cutoff (+60 eV): {calc_cutoff} eV"
    else:
        calc_cutoff = np.ceil(max_base_cutoff / 10.0) * 10
        cutoff_log_msg = f"Exact Cutoff (No Buffer): {calc_cutoff} eV"

    for sym in symbols:
        source = OTFG_STRINGS.get(sym, "C19")
        logger.info(f"  [{sym}] Injected OTFG String: {source}")

    logger.info(f"  -> Target Quality: {cutoff_quality.upper()} | Highest Base: {max_base_cutoff:.1f} eV")
    logger.info(f"  -> K-Point Spacing: {kpt_spacing} 1/A")
    logger.info(f"  -> {cutoff_log_msg}")

    with open(cell_path, "r", encoding="utf-8") as f: lines = f.read().splitlines()
    
    clean_lines, in_b, cur_b = [], False, ""
    for l in lines:
        l_low = l.strip().lower()
        if "fix_com" in l_low or "kpoints_mp_grid" in l_low or "kpoints_mp_spacing" in l_low: continue
        
        if not in_b:
            for b in ["kpoints_list", "ionic_constraints", "external_efield", "external_pressure", "species_mass", "species_pot"]:
                if l_low.startswith(f"%block {b}"): in_b, cur_b = True, b; break
            if in_b: continue
        elif l_low.startswith(f"%endblock {cur_b}"): in_b = False; continue
        if not in_b: clean_lines.append(l)

    cell_text, total_spin = inject_quantum_spins("\n".join(clean_lines) + "\n")
    
    k_block = [f"KPOINTS_MP_SPACING {kpt_spacing}\n"]
    cell_text = cell_text.replace("%ENDBLOCK POSITIONS_FRAC\n", "%ENDBLOCK POSITIONS_FRAC\n\n" + "\n".join(k_block))

    tail = ["FIX_COM : false", "", "%BLOCK EXTERNAL_PRESSURE", "    0.0  0.0  0.0\n         0.0  0.0\n              0.0", "%ENDBLOCK EXTERNAL_PRESSURE", ""]
        
    tail.extend(["%BLOCK SPECIES_MASS"] + [f"      {s}    {atomic_masses[atomic_numbers.get(s)]:.10f}" for s in symbols] + ["%ENDBLOCK SPECIES_MASS", ""])
    
    tail.extend(["%BLOCK SPECIES_POT"] + [f"      {s}  {OTFG_STRINGS.get(s, 'C19')}" for s in symbols] + ["%ENDBLOCK SPECIES_POT", ""])
    
    with open(cell_path, "w", encoding="utf-8") as f: f.write(cell_text + "\n".join(tail))

    params = config.get("castep_base_params", {}).copy()
    
    params.update({
        "cut_off_energy": f"{calc_cutoff:.15f}",
        "nextra_bands": calc_nextra, 
        "geom_modulus_est": f"{props['modulus_est']} GPa"
    })
    
    if props["charge"] != 0: params["charge"] = props["charge"]
    if total_spin > 0: params.update({"spin_polarized": True, "spin": int(total_spin)})
    
    if props["complexity"] > 2.0:
        params["geom_max_iter"] = 1500
        params["max_scf_cycles"] = 900

    params["fix_occupancy"] = False
    params["metals_method"] = "dm"
    params["smearing_width"] = "0.010000000000000"

    with open(os.path.join(work_dir, f"{seed}.param"), "w") as f:
        for k, v in params.items(): f.write(f"{k} : {'true' if v is True else 'false' if v is False else v}\n")

    logger.info(f"OK: {seed} | Formula: {props['formula']} (Z={props['z_factor']}) | Explicit Cutoff: {calc_cutoff} eV")


if __name__ == "__main__":
    if getattr(sys, 'frozen', False):
        BASE_DIR = os.path.dirname(sys.executable)
    else:
        BASE_DIR = os.path.dirname(os.path.abspath(__file__))
        
    default_config = os.path.join(BASE_DIR, "config.json")

    parser = argparse.ArgumentParser(description="Enterprise CASTEP Pipeline")
    parser.add_argument("--config", default=default_config, help="Path to config file")
    parser.add_argument("--cif_dir", type=str, help="Override path to input CIF directory")
    parser.add_argument("--precision", type=str.lower, choices=['coarse', 'medium', 'fine', 'ultrafine'], default='ultrafine', help="Set the CASTEP cutoff energy precision (COARSE, MEDIUM, FINE, or ULTRAFINE). Defaults to ULTRAFINE.")
    args = parser.parse_args()

    if not os.path.exists(args.config):
        logger.error(f"CRITICAL: Cannot find config.json at {args.config}")
        sys.exit(1)

    config = load_config(args.config)

    if args.cif_dir:
        if "directories" not in config:
            config["directories"] = {}
        config["directories"]["input_cifs"] = args.cif_dir

    if args.precision:
        if "heuristics" not in config:
            config["heuristics"] = {}
        config["heuristics"]["cutoff_quality"] = args.precision

    cifs = sorted(glob.glob(os.path.join(os.path.expanduser(config["directories"]["input_cifs"]), "*.cif")))

    if not cifs:
        logger.warning(f"No CIF files found in {config['directories']['input_cifs']}")

    for cif in cifs:
        try:
            process_one_cif(cif, config)
        except Exception as e:
            logger.error(f"FAIL: {os.path.basename(cif)} -> {e}")

