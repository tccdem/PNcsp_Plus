import os
import pandas as pd
from ase.io import read, write
import argparse

def general_relaxer(atoms, calculator, logfile, fmax=0.05, steps=500, CellFilter=True):
    from ase.optimize.fire import FIRE
    from ase.constraints import ExpCellFilter
    print("Relaxer starts...")

    atoms.calc = calculator

    
    if(CellFilter==True):
        atoms = ExpCellFilter(atoms)
        dyn = FIRE(atoms,logfile=logfile)
        dyn.run(fmax=fmax, steps=steps)
        if not dyn.converged():
            with open(logfile, "a") as f:
                f.write("Warning: Not Converged!\n")
        else:
            with open(logfile, "a") as f:
                f.write("Converged!\n")
        return atoms.atoms,  atoms.atoms.get_potential_energy()
    else:
        dyn = FIRE(atoms,logfile=logfile)
        dyn.run(fmax=fmax, steps=steps)
        # return ase_to_atoms(ase_atoms.atoms), ase_atoms.get_potential_energy()
        if not dyn.converged():
            with open(logfile, "a") as f:
                f.write("Warning: Not Converged!\n")
                print("Warning: Not Converged!")
        else:
            with open(logfile, "a") as f:
                f.write("Converged!\n")
                print("Converged!")
        return atoms, atoms.get_potential_energy()

def MACE_calc(formula,path):
    from mace.calculators import mace_mp
    import os
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.io.cif import CifWriter



    path_results=os.path.join(path,"output_"+formula+"/Calc_report/MACE/CheckNew_report/Best_Structures")
    dest_path=os.path.join(path_results,"opt")

    if not os.path.exists(dest_path):
        os.makedirs(dest_path)
    if not os.path.exists(path_results):
        exit(0)

    logfile=os.path.join(dest_path,"RelaxationReport.txt")
    if os.path.exists(logfile):
        os.remove(logfile)

    struct_list_init=[]
    struc_name_list=os.listdir(path_results)
    if(struc_name_list==[]):
        exit(0)

    print("\nMACE is activated...\n")
    calc = mace_mp(
        model="medium-mpa-0",        # “small” | “medium” | “large” | medium-mpa-0 etc
        dispersion=False,      # no empirical dispersion correction
        default_dtype="float64",
        device="cpu"          # or "cpu"
    )

    with open(logfile, "a") as f:
        for struc_name in struc_name_list:

            if("opt" in struc_name):
                continue

            print(f"{struc_name}")
            f.write(struc_name+"\n"+"*****************\n")
            f.flush()
            target_path=os.path.join(path_results,struc_name)
            struct = read(target_path)
            atom_num=struct.get_global_number_of_atoms()
            final_structure, final_energy = general_relaxer(atoms=struct, calculator=calc,logfile=logfile,CellFilter=True)

            print(f"Predicted Energy= {final_energy/atom_num:.4f} eV/atom")

            # write(os.path.join(path_results,struc_name.replace(".cif","_opt.cif")), final_structure)
            pmg_structure = AseAtomsAdaptor.get_structure(final_structure)
            writer = CifWriter(pmg_structure, symprec=1e-2, angle_tolerance=5, refine_struct=True)
            
            writer.write_file(os.path.join(dest_path,struc_name.replace(".cif","_opt.cif")))
    print("Terminated Sucessfully!")  

def main():
    parser = argparse.ArgumentParser(prog="Relaxer",description= "Relaxes with MACE")
    parser.add_argument('formula')
    parser.add_argument('-out','--output_dir',default=".",help="Sets output directory. Enter full path. (default: current directory).")
 
    args = parser.parse_args()

    formula=args.formula

    if(args.output_dir[-1]=="/"):
        data_path=args.output_dir[:-1]
    else:
        data_path=args.output_dir

    MACE_calc(formula,data_path)

if __name__=='__main__':
    main()