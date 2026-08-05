"""
PNcsp+ : Interpretable Crystal Structure Prediction Framework
"""

# 1. Package Metadata
__version__ = "0.1.0"
__author__ = "Cem Oran"

# 2. Expose core functions/classes at the root package level
from .PNcsp import run, get_Neig, get_Symbol, get_PN, convert_formula
from PNcsp_Plus.db import DBsearch

# 3. Control public exports for wildcard imports (`from PNcsp_Plus import *`)
__all__ = [
    "run",
    "get_Neig",
    "get_Symbol",
    "get_PN",
    "convert_formula",
    "__version__",
]