![Graphical_Abstract_3](https://github.com/user-attachments/assets/cf5be4d4-2d3e-4ccc-9aea-df1c756dafca
)
[![DOI](https://img.shields.io/badge/DOI-10.1021/acs.jctc.6c00044-blue?style=flat-square)](https://pubs.acs.org/doi/full/10.1021/acs.jctc.6c00044)
# PNcsp+

This repository contains the AI assited version of the code implementation of the CSP method, as discussed in the paper [PNcsp+: A Periodic Number-Based Crystal Structure Prediction Method Enhanced by Machine Learning](https://pubs.acs.org/doi/full/10.1021/acs.jctc.6c00044), developed by the [TCCDEM](https://github.com/tccdem/) team. Currently MACE, M3GNet and ALIGNN-FF are suported.

This work introduces a novel strategy for crystal structure prediction founded upon the principle of chemical similarity.  Our method uses Mendeleev's Periodic Number (PN) as a quantitative measure of substitutability to identify potential crystal structures for unexplored chemical systems. Representation of the workflow for predicting stable materials based on PN similarity is shown below. 

![Workflow](https://github.com/user-attachments/assets/37c070da-f9c4-4f12-8b44-00a3f08d9837)

PNcsp is a generaliezed version of this workflow. It scans the [OQMD](https://www.oqmd.org/) for a desired chemical system for given order of neighbors and proposes initial crystal structures that are ready for further analysis.

**PLEASE NOTE:** the PNcsp code is under active development. Bug reports are welcomed in the GitHub issues!

### Installation & Usage
This program is based on Python 3 under Anaconda. 

1) Clone the repository.
2) Open terminal and locate the directory.
3) Install reuqirements or run:
```bash
   pip install -r requirements.txt
```
4) Help page:
```bash
   python PNcsp.py -h
```
5) Run the Python code:
```bash
  python PNcsp.py <formula> -n <neighbor_order>  -f <energy_filter> -c <calculator>
```

- `-n`, `--neighbor`  
  Order of neighbors to be considered in the similarity search.

- `-f`, `--filter`  
  Selected neighbors are limited to those below the energy filter value.  
  **Default:** `0` (unit: eV/atom).  
  Use `"none"` to disable filtering.

  - `-db`, `--database`  
  Sets the data source: `OQMD`, `MP`, or `MPDS`.  
  **Default:** `OQMD`.

- `-calc`, `--calculator`  
  Selects calculator:  
  `M3GNet`, `ALIGNN`, `MACE`, or `Ensemble`.  
  **Default:** `None`.

- `-out`, `--output_dir`  
  Sets output directory. Provide the full path.  
  **Default:** current directory.

- `-o`, `--online`  
  Enables online (`True`) or offline (`False`) search in OQMD.  
  **Default:** `False`.

  For offline search, download and set up the offline OQMD database:  
  https://oqmd.org/download/

- `-t`, `--time_sleep`  
  Sets sleep time between queries for online search in OQMD.  
  Excessive number of queries may cause the server to halt.  
  **Default:** `"none"`.

- `--BlockSearch`  
  Blocks the search stage.  
  Use this flag if you want to run only the calculator without the search feature.

- `--Relax`  
  Performs structure relaxation before ML evaluation.




Created prototypes are shown in "output" folder in current directory.

### Example usage
```bash
python Similarity.py Na2Cl1 -n 3 -f 0.1 -c MACE -db ./Outputs
```
