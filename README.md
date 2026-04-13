![Graphical_Abstract_3](https://github.com/user-attachments/assets/cf5be4d4-2d3e-4ccc-9aea-df1c756dafca
)
[![DOI](https://img.shields.io/badge/DOI-10.1021/acs.jctc.6c00044-blue?style=flat-square)](https://pubs.acs.org/doi/full/10.1021/acs.jctc.6c00044)
# PNcsp+

This repository contains the AI assited version of the code implementation of the CSP method, as discussed in the paper [PNcsp+: A Periodic Number-Based Crystal Structure Prediction Method Enhanced by Machine Learning](https://pubs.acs.org/doi/full/10.1021/acs.jctc.6c00044), developed by the [TCCDEM](https://github.com/tccdem/) team. Currently MACE, M3GNet and ALIGNN-FF are suported.

This work introduces a novel strategy for crystal structure prediction founded upon the principle of chemical similarity.  Our method uses Mendeleev's Periodic Number (PN) as a quantitative measure of substitutability to identify potential crystal structures for unexplored chemical systems. Representation of the workflow for predicting stable materials based on PN similarity is shown below. 

![Workflow](https://github.com/user-attachments/assets/37c070da-f9c4-4f12-8b44-00a3f08d9837)

PNcsp+ with its default configuration scans a local instance of the [OQMD](https://www.oqmd.org/) for a desired chemical system for given order of neighbors and proposes crystal structures that are ready for further analysis. To set up OQMD locally, please refer to the [qmpy documantation](https://static.oqmd.org/static/docs/getting_started.html#setting-up-the-database).

**PLEASE NOTE:** the PNcsp code is under active development. Bug reports are welcomed in the GitHub issues!

### Installation & Usage
This program is based on Python 3 under Anaconda. 

1) Clone the repository.
2) Open terminal and locate the directory.
3) Install reuqirements (check also all_requirements.txt):
```bash
   pip install -r requirements.txt
```
4) Help page:
```bash
   python PNcsp.py -h
```
5) Run the Python code (minimal):
```bash
  python PNcsp.py <formula> -n <neighbor_order>  -f <energy_filter> -c <calculator>
```

- `-n`, `--neighbor`  
  Order of neighbors to be considered in the similarity search.
  **Default:** `1`.

- `-f`, `--filter`  
  Selected neighbors are limited to those below the energy filter value.  
  **Default:** `0` (unit: eV/atom).  
  Use `"none"` to disable filtering.

  - `-db`, `--database`  
  Sets the data source: `OQMD`, `MP`, or `MPDS`.  
  **Default:** `OQMD`.

- `-calc`, `--calculator`  
  Selects calculator:  
  `M3GNet`, `ALIGNN`, `MACE`, or `ensemble`.  
  **Default:** `None`.

- `-out`, `--output_dir`  
  Sets output directory. Provide the full path.  
  **Default:** current directory.

- `-on`, `--online`  
  Enables online (`True`) or offline (`False`) search in OQMD.  
  **Default:** `False`.

  For offline search, download and set up the offline OQMD database:  
  https://oqmd.org/download/

- `-ts`, `--time_sleep`  
  Sets sleep time between queries for online search in OQMD.  
  Excessive number of queries may cause the server to halt.  
  **Default:** `"none"`.

- `--BlockSearch`  
  Blocks the search stage.  
  Use this flag if you want to run only the calculator without the search feature.

- `--Relax`  
  Performs structure relaxation before ML evaluation.
- `--CheckNew`
   Check if found structures have been already reported in OQMD and MP.
- `-top_n`, `--top_n_new`
  Copies the top-n evaluated new structures, ranked by the GNN evaluation, to the Best_Structures folder if available [int, "all", "none"]. (default: none).
  Use this option together with --CheckNew, or in a subsequent run after a run performed with --CheckNew.
- `-top_c`, `--top_n_calc`
  Copies the top-n evaluated structures, ranked by the GNN evaluation, to the Best_Structures folder if available [int, "all", "none"].(default: none). 
  Use this option together with --CheckNew, or in a subsequent run after a run performed with --CheckNew.




Created prototypes are shown in "output" folder in current directory.

### Example usage
```bash
python PNcsp.py Na2Cl1 -n 3 -f 0.1 -c MACE -out ./output_dir
```

