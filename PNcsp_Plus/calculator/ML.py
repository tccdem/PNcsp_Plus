import os
import pandas as pd
from ase.io import read, write

def general_relaxer(atoms, calculator, logfile, fmax=0.05, steps=500, CellFilter=True):
    from ase.optimize.fire import FIRE
    from ase.constraints import ExpCellFilter
    print("Relaxer starts...")
    # print(f"Number of atoms in CIF: {atoms.num_atoms}")
    # print(f"Number of atoms in CIF: {atoms.get_number_of_atoms()}")
    # ase_atoms = atoms.ase_converter()
    # ase_atoms = atoms
    atoms.calc = calculator
    # if not relax:
    #     return ase_atoms, ase_atoms.get_potential_energy()
    
    if(CellFilter==True):
        atoms = ExpCellFilter(atoms)
        dyn = FIRE(atoms,logfile=logfile)
        dyn.run(fmax=fmax, steps=steps)
        # return ase_to_atoms(ase_atoms.atoms), ase_atoms.get_potential_energy()
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

def MACE_calc(formula,path,path_results,relax):
    from mace.calculators import mace_mp
    # from ase.optimize.fire import FIRE
    # from ase.constraints import ExpCellFilter
    # import ase
    

    # atoms = read("./Si2_modified.cif")
    # atom_num=atoms.get_global_number_of_atoms()

    print("\nMACE is activated...\n")
    if not os.path.exists(path_results):
        os.makedirs(path_results)
    # comp_list=os.listdir(path)
    # for comp in comp_list:
    #     if not "output_" in comp:
    #         continue
    #     print(comp)
    struct_list_init=[]
    neigh_list=os.listdir(path)
    if(neigh_list==[]):
        exit(0)
    for neigh in neigh_list:
        if "Neigh" not in neigh:
            continue
        # if neigh!="1_Neigh":
        #     continue
        sym_list=os.listdir(path+neigh)
        if(sym_list==[]):
            continue
        for sym in sym_list:
            proto_list=os.listdir(path+neigh+"/"+sym)
            if(proto_list==[]):
                continue
            for proto in proto_list:
                target_path=path+neigh+"/"+sym+"/"+proto
                # print(comp,neigh,sym,proto)
                # struct=Structure.from_file(target_path)
                # struct = Atoms.from_cif(target_path)
                struct = read(target_path)
                struct_list_init.append([struct,formula,neigh,sym,proto])
                # print(proto)
    print("\nTotal number of systems:",len(struct_list_init))
    print("------------------------")
    
    
    final_structure_list=[]


    # def general_relaxer(atoms, calculator, fmax=0.05, steps=500, CellFilter=True):
    #     print("Relaxer starts...")
    #     # print(f"Number of atoms in CIF: {atoms.num_atoms}")
    #     # print(f"Number of atoms in CIF: {atoms.get_number_of_atoms()}")
    #     # ase_atoms = atoms.ase_converter()
    #     # ase_atoms = atoms
    #     atoms.calc = calculator
    #     # if not relax:
    #     #     return ase_atoms, ase_atoms.get_potential_energy()
        
    #     if(CellFilter==True):
    #         atoms = ExpCellFilter(atoms)
    #         dyn = FIRE(atoms,logfile=logfile)
    #         dyn.run(fmax=fmax, steps=steps)
    #         # return ase_to_atoms(ase_atoms.atoms), ase_atoms.get_potential_energy()
    #         if not dyn.converged():
    #             with open(logfile, "a") as f:
    #                 f.write("Warning: Not Converged!\n")
    #         else:
    #             with open(logfile, "a") as f:
    #                 f.write("Converged!\n")
    #         return atoms.atoms,  atoms.atoms.get_potential_energy()
    #     else:
    #         dyn = FIRE(atoms,logfile=logfile)
    #         dyn.run(fmax=fmax, steps=steps)
    #         # return ase_to_atoms(ase_atoms.atoms), ase_atoms.get_potential_energy()
    #         if not dyn.converged():
    #             with open(logfile, "a") as f:
    #                 f.write("Warning: Not Converged!\n")
    #                 print("Warning: Not Converged!")
    #         else:
    #             with open(logfile, "a") as f:
    #                 f.write("Converged!\n")
    #                 print("Converged!")
    #         return atoms, atoms.get_potential_energy()

    calc = mace_mp(
        model="medium-mpa-0",        # “small” | “medium” | “large” | medium-mpa-0 etc
        dispersion=False,      # no empirical dispersion correction
        default_dtype="float64",
        device="cpu"          # or "cpu"
    )

    if relax==True:
        logfile=path_results+"RelaxationReport.txt"
        if os.path.exists(logfile):
            os.remove(logfile)
        for struct in struct_list_init:
            print("\n"+struct[-1])
            # if not "Pm-3m" in struct[-1]:
            #     continue
            with open(logfile, "a") as f:
                f.write(struct[-1]+"\n"+"*****************\n")
            # print("Initial structure:\n", struct[0])
            # atoms = general_relaxer(atoms=struct, calculator=calc)
            # atom_num=struct[0].get_number_of_atoms()
            atom_num=struct[0].get_global_number_of_atoms()
            final_structure, final_energy = general_relaxer(atoms=struct[0], calculator=calc,logfile=logfile,CellFilter=True)
            # print("Initial structure:\n", struct[0])
            # print("Relaxed structure:\n", final_structure)
            print(f"Predicted Energy= {final_energy/atom_num:.4f} eV/atom")
            final_structure_list.append({"CIF_Name":struct[4].replace(".cif",""),"Sym":struct[3],"Neigh_Order":struct[2],"Energy":final_energy/atom_num,"Target":struct[1],"Struc_Init":struct[0],"Struc_Final":final_structure})

        
            if not os.path.exists(path_results+struct[4]+"_opt.cif"):
                # final_structure.to_file(path_results+struct[4]+"_opt.cif")
                write(path_results+struct[4]+"_MACE.cif", final_structure)
    
        df=pd.DataFrame(final_structure_list)
        # df['Energy'] = pd.to_numeric(df['Energy'], errors='coerce').round(3)
        df['Energy'] = df['Energy'].round(3)
        df=df.sort_values(by=["Energy"])
        df[["CIF_Name","Neigh_Order","Energy"]].to_csv(path_results+"MACE_"+final_structure_list[0]["Target"]+'_relaxed_all.csv', index=False)

        df=df.drop_duplicates(subset=["Sym"])
        df[["CIF_Name","Neigh_Order","Energy"]].to_csv(path_results+"MACE_"+final_structure_list[0]["Target"]+'_relaxed_best.csv', index=False)


    else:
        for struct in struct_list_init:
            print("\n"+struct[-1])
            # print("Initial structure:\n", struct[0])
            # atoms = general_relaxer(atoms=struct, calculator=calc)
            # atom_num=struct[0].get_number_of_atoms()
            atom_num=struct[0].get_global_number_of_atoms()

            struct[0].calc = calc
            final_energy = struct[0].get_potential_energy()
                    
            # print("Initial structure:\n", struct[0])
            # print("Relaxed structure:\n", final_structure)
            print(f"Predicted Energy= {final_energy/atom_num:.4f} eV/atom")
            final_structure_list.append({"CIF_Name":struct[4].replace(".cif",""),"Sym":struct[3],"Neigh_Order":struct[2],"Energy":final_energy/atom_num,"Target":struct[1],"Struc_Init":struct[0]})
    
        df=pd.DataFrame(final_structure_list)
        # df['Energy'] = pd.to_numeric(df['Energy'], errors='coerce').round(3)
        df['Energy'] = df['Energy'].round(3)
        df=df.sort_values(by=["Energy"])
        df[["CIF_Name","Neigh_Order","Energy"]].to_csv(path_results+"MACE_"+final_structure_list[0]["Target"]+'_all.csv', index=False)

        df=df.drop_duplicates(subset=["Sym"])
        df[["CIF_Name","Neigh_Order","Energy"]].to_csv(path_results+"MACE_"+final_structure_list[0]["Target"]+'_best.csv', index=False)

    print("Done")
    print("------------------------")
    print("Terminated Sucessfully!")  
    
    # # atoms.calc = calc
    # final_structure, final_energy = general_relaxer(atoms=atoms,calculator=calc,relax=True,CellFilter=False)

    # # ────────────────────────────────────────────────────────────────────────────────
    # # 3) Compute & report energy, then save CIF
    # # ────────────────────────────────────────────────────────────────────────────────
    # energy = atoms.get_potential_energy()
    # # print(f"Si2 potential energy: {energy:.6f} eV per atom")
    # print(f"Initial Predicted Energy= {energy/atom_num:.4f} eV/atom")
    # print(f"Final Predicted Energy= {final_energy/atom_num:.4f} eV/atom")
    # write("./Si2_opt.cif", final_structure)

def ALIGNN_calc(formula,path,path_results,relax):
    from alignn.ff.ff import AlignnAtomwiseCalculator, default_path
    from jarvis.core.atoms import Atoms, ase_to_atoms
    # from ase.constraints import ExpCellFilter
    # from ase.optimize.fire import FIRE
    # import matplotlib.pyplot as plt
    # import numpy as np

    from pymatgen.core import Lattice, Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    print("\nALIGNN-FF is activated...\n")
    if not os.path.exists(path_results):
        os.makedirs(path_results)
    # comp_list=os.listdir(path)
    # for comp in comp_list:
    #     if not "output_" in comp:
    #         continue
    #     print(comp)
    struct_list_init=[]
    neigh_list=os.listdir(path)
    if(neigh_list==[]):
        exit(0)
    for neigh in neigh_list:
        if "Neigh" not in neigh:
            continue
        # if neigh!="1_Neigh":
        #     continue
        sym_list=os.listdir(path+neigh)
        if(sym_list==[]):
            continue
        for sym in sym_list:
            proto_list=os.listdir(path+neigh+"/"+sym)
            if(proto_list==[]):
                continue
            for proto in proto_list:
                target_path=path+neigh+"/"+sym+"/"+proto
                # print(comp,neigh,sym,proto)
                # struct=Structure.from_file(target_path)
                # struct = Atoms.from_cif(target_path)
                struct = read(target_path)
                struct_list_init.append([struct,formula,neigh,sym,proto])
                # print(proto)
    print("\nTotal number of systems:",len(struct_list_init))
    print("------------------------")
    
    
    final_structure_list=[]

    # # Relaxer function
    # def general_relaxer(atoms, calculator, fmax=0.05, steps=500, CellFilter=True):
    #     print("Relaxer starts...")
    #     # print(f"Number of atoms in CIF: {atoms.num_atoms}")
    #     # print(f"Number of atoms in CIF: {atoms.get_number_of_atoms()}")
    #     # ase_atoms = atoms.ase_converter()
    #     # ase_atoms = atoms
    #     atoms.calc = calculator
    #     # if not relax:
    #     #     return ase_atoms, ase_atoms.get_potential_energy()
        
    #     if(CellFilter==True):
    #         atoms = ExpCellFilter(atoms)
    #         dyn = FIRE(atoms,logfile=logfile)
    #         dyn.run(fmax=fmax, steps=steps)
    #         # return ase_to_atoms(ase_atoms.atoms), ase_atoms.get_potential_energy()
    #         if not dyn.converged():
    #             with open(logfile, "a") as f:
    #                 f.write("Warning: Not Converged!\n")
    #         else:
    #             with open(logfile, "a") as f:
    #                 f.write("Converged!\n")
    #         return atoms.atoms,  atoms.atoms.get_potential_energy()
    #     else:
    #         dyn = FIRE(atoms,logfile=logfile)
    #         dyn.run(fmax=fmax, steps=steps)
    #         # return ase_to_atoms(ase_atoms.atoms), ase_atoms.get_potential_energy()
    #         if not dyn.converged():
    #             with open(logfile, "a") as f:
    #                 f.write("Warning: Not Converged!\n")
    #                 print("Warning: Not Converged!")
    #         else:
    #             with open(logfile, "a") as f:
    #                 f.write("Converged!\n")
    #         return atoms, atoms.get_potential_energy()

    
    # ALIGNN Calculator
    model_path = default_path()
    # model_path = "/home/cem/alignn/alignn/ff/jv_formation_energy_peratom_alignn"
    calc = AlignnAtomwiseCalculator(path=model_path)

    # # Load CIF
    # atoms = Atoms.from_cif("./Si2.cif")

    # Main logic
    if relax==True:
        logfile=path_results+"RelaxationReport.txt"
        if os.path.exists(logfile):
            os.remove(logfile)
        for struct in struct_list_init:
            print("\n"+struct[-1])
            with open(logfile, "a") as f:
                f.write(struct[-1]+"\n"+"*****************\n")
            # print("Initial structure:\n", struct[0])
            # atoms = general_relaxer(atoms=struct, calculator=calc)
            # atom_num=struct[0].get_number_of_atoms()
            atom_num=struct[0].get_global_number_of_atoms()
            final_structure, final_energy = general_relaxer(atoms=struct[0], calculator=calc, logfile=logfile, CellFilter=True)
            # print("Initial structure:\n", struct[0])
            # print("Relaxed structure:\n", final_structure)
            print(f"Predicted Energy= {final_energy/atom_num:.4f} eV/atom")
            final_structure_list.append({"CIF_Name":struct[4].replace(".cif",""),"Sym":struct[3],"Neigh_Order":struct[2],"Energy":final_energy/atom_num,"Target":struct[1],"Struc_Init":struct[0],"Struc_Final":final_structure})

            
            if not os.path.exists(path_results+struct[4]+"_opt.cif"):
                # final_structure.to_file(path_results+struct[4]+"_opt.cif")
                write(path_results+struct[4]+"_ALIGNN.cif", final_structure)
        
        df=pd.DataFrame(final_structure_list)
        # df['Energy'] = pd.to_numeric(df['Energy'], errors='coerce').round(3)
        df['Energy'] = df['Energy'].round(3)
        df=df.sort_values(by=["Energy"])
        df[["CIF_Name","Neigh_Order","Energy"]].to_csv(path_results+"ALIGNNFF_"+final_structure_list[0]["Target"]+'_relaxed_all.csv', index=False)

        df=df.drop_duplicates(subset=["Sym"])
        df[["CIF_Name","Neigh_Order","Energy"]].to_csv(path_results+"ALIGNNFF_"+final_structure_list[0]["Target"]+'_relaxed_best.csv', index=False)

    else:
        for struct in struct_list_init:
            print("\n"+struct[-1])
            # print("Initial structure:\n", struct[0])
            # atoms = general_relaxer(atoms=struct, calculator=calc)
            # atom_num=struct[0].get_number_of_atoms()
            atom_num=struct[0].get_global_number_of_atoms()

            struct[0].calc = calc
            final_energy = struct[0].get_potential_energy()
                    
            # print("Initial structure:\n", struct[0])
            # print("Relaxed structure:\n", final_structure)
            print(f"Predicted Energy= {final_energy/atom_num:.4f} eV/atom")
            final_structure_list.append({"CIF_Name":struct[4].replace(".cif",""),"Sym":struct[3],"Neigh_Order":struct[2],"Energy":final_energy/atom_num,"Target":struct[1],"Struc_Init":struct[0]})

        
        df=pd.DataFrame(final_structure_list)
        # df['Energy'] = pd.to_numeric(df['Energy'], errors='coerce').round(3)
        df['Energy'] = df['Energy'].round(3)
        df=df.sort_values(by=["Energy"])
        df[["CIF_Name","Neigh_Order","Energy"]].to_csv(path_results+"ALIGNNFF_"+final_structure_list[0]["Target"]+'_all.csv', index=False)

        df=df.drop_duplicates(subset=["Sym"])
        df[["CIF_Name","Neigh_Order","Energy"]].to_csv(path_results+"ALIGNNFF_"+final_structure_list[0]["Target"]+'_best.csv', index=False)
    print("Done")
    print("------------------------")
    print("Terminated Sucessfully!")

def M3GNet_calc(formula,path,path_results,relax):
    import warnings

    from m3gnet.models import Relaxer
    from pymatgen.core import Lattice, Structure
    
    from m3gnet.models import M3GNet

    # for category in (UserWarning, DeprecationWarning):
    #     warnings.filterwarnings("ignore", category=category, module="tensorflow")

    # import tensorflow as tf
    # tf.get_logger().setLevel('ERROR')


    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    # path="./Unknown_Systems_Data/Negative/"
    
    print("\nM3GNET is activated...\n")
    if not os.path.exists(path_results):
        os.makedirs(path_results)
    # comp_list=os.listdir(path)
    # for comp in comp_list:
    #     if not "output_" in comp:
    #         continue
    #     print(comp)
    struct_list_init=[]
    neigh_list=os.listdir(path)
    if(neigh_list==[]):
        exit(0)
    for neigh in neigh_list:
        if "Neigh" not in neigh:
            continue
        # if neigh!="1_Neigh":
        #     continue
        sym_list=os.listdir(path+neigh)
        if(sym_list==[]):
            continue
        for sym in sym_list:
            proto_list=os.listdir(path+neigh+"/"+sym)
            if(proto_list==[]):
                continue
            for proto in proto_list:
                target_path=path+neigh+"/"+sym+"/"+proto
                # print(comp,neigh,sym,proto)
                struct=read(target_path)
                struct_list_init.append([struct,formula,neigh,sym,proto])
                # print(proto)
    print("\nTotal number of systems:",len(struct_list_init))
    print("------------------------")

    if(relax==True):
        print("Relaxer starts...")

        logfile=path_results+"RelaxationReport.txt"
        if os.path.exists(logfile):
            os.remove(logfile)
        
        final_structure_list=[]
        relaxer=Relaxer(relax_cell=True)
        for struct in struct_list_init:
            print("\n"+struct[-1])
            with open(logfile, "a") as f:
                f.write(struct[-1]+"\n"+"*****************\n")

            relax_results = relaxer.relax(struct[0], fmax=0.05, steps=500, verbose=True,logfile=logfile)

            # if not relaxer.opt_class.converged:
            #     with open(logfile, "a") as f:
            #         f.write("Warning: Not Converged!\n")
            # else:
            #     with open(logfile, "a") as f:
            #         f.write("Converged!\n")

            atom_num=struct[0].get_global_number_of_atoms()
            final_structure=relax_results['final_structure']
            final_energy_per_atom=float(relax_results['trajectory'].energies[-1] / atom_num)
            print(f"Predicted Energy= {final_energy_per_atom:.4f} eV/atom")
            final_structure_list.append({"CIF_Name":struct[4].replace(".cif",""),"Sym":struct[3],"Neigh_Order":struct[2],"Energy":final_energy_per_atom,"Target":struct[1],"Struc_Init":struct[0],"Struc_Final":final_structure,})
            if not os.path.exists(path_results+struct[4]+"_opt.cif"):
                final_structure.to_file(path_results+struct[4]+"_M3Gnet.cif")
            # print(f"Relaxed lattice parameter is {final_structure.lattice.abc[0]:.3f} Å")
            # print(f"Final energy is {final_energy_per_atom:.3f} eV/atom")
        df=pd.DataFrame(final_structure_list)
        df['Energy'] = df['Energy'].round(3)
        df=df.sort_values(by=["Energy"])
        # df[["CIF_Name","Sym","Neigh_Order","Energy","Target"]].to_csv(path_results+"M3Gnet_"+final_structure_list[0]["Target"]+'_relaxed_all.csv', index=False)
        df[["CIF_Name","Neigh_Order","Energy"]].to_csv(path_results+"M3Gnet_"+final_structure_list[0]["Target"]+'_relaxed_all.csv', index=False)

        df=df.drop_duplicates(subset=["Sym"])
        # df[["CIF_Name","Sym","Neigh_Order","Energy","Target"]].to_csv(path_results+"M3Gnet_"+final_structure_list[0]["Target"]+'_relaxed_best.csv', index=False)
        df[["CIF_Name","Neigh_Order","Energy"]].to_csv(path_results+"M3Gnet_"+final_structure_list[0]["Target"]+'_relaxed_best.csv', index=False)
    else:
        final_structure_list=[]
        model = M3GNet.load()  # defaults to "MP-2021.2.8-EFS" :contentReference
        for struct in struct_list_init:
            print("\n"+struct[-1])
            final_structure=struct[0]

            atom_num=struct[0].get_global_number_of_atoms()
            final_energy_per_atom = float(model.predict_structure(struct[0])/atom_num)
            print(f"Predicted Energy= {final_energy_per_atom:.4f} eV/atom")
            final_structure_list.append({"CIF_Name":struct[4].replace(".cif",""),"Sym":struct[3],"Neigh_Order":struct[2],"Energy":final_energy_per_atom,"Target":struct[1],"Struc_Init":struct[0]})
            # if not os.path.exists(path_results+struct[4]+"_opt.cif"):
            #     final_structure.to_file(path_results+struct[4]+"_opt.cif")
            # print(f"Relaxed lattice parameter is {final_structure.lattice.abc[0]:.3f} Å")
            # print(f"Final energy is {final_energy_per_atom:.3f} eV/atom")
            df=pd.DataFrame(final_structure_list)
            df['Energy'] = df['Energy'].round(3)
            df=df.sort_values(by=["Energy"])
            # df[["CIF_Name","Sym","Neigh_Order","Energy","Target"]].to_csv(path_results+"M3Gnet_"+final_structure_list[0]["Target"]+'_all.csv', index=False)
            df[["CIF_Name","Neigh_Order","Energy"]].to_csv(path_results+"M3Gnet_"+final_structure_list[0]["Target"]+'_all.csv', index=False)

            df=df.drop_duplicates(subset=["Sym"])
            # df[["CIF_Name","Sym","Neigh_Order","Energy","Target"]].to_csv(path_results+"M3Gnet_"+final_structure_list[0]["Target"]+'_best.csv', index=False)
            df[["CIF_Name","Neigh_Order","Energy"]].to_csv(path_results+"M3Gnet_"+final_structure_list[0]["Target"]+'_best.csv', index=False)
    print("Done")
    print("------------------------")
    print("Terminated Sucessfully!")

def majority_vote(formula,path_alignn,path_mace,path_m3gnet,path_results, N_model=3):
    if not os.path.exists(path_results):
        os.makedirs(path_results)

    if N_model==3:
        df_alignn=pd.read_csv(path_alignn+"ALIGNNFF_"+formula+"_all.csv")
        df_m3gnet=pd.read_csv(path_m3gnet+"M3Gnet_"+formula+"_all.csv")
        df_mace=pd.read_csv(path_mace+"MACE_"+formula+"_all.csv")

        # merged_df=pd.merge(df_alignn, df_m3gnet, on='CIF_Name', how='inner')
        merged_df = pd.merge(df_alignn, df_m3gnet, on=['CIF_Name', 'Neigh_Order'], suffixes=('_alignn', '_m3gnet'))

        merged_df = pd.merge(merged_df, df_mace,on=['CIF_Name', 'Neigh_Order'])
        merged_df = merged_df.rename(columns={'Energy':'Energy_mace'})

        merged_df['Mean_Energy'] = merged_df[['Energy_alignn', 'Energy_m3gnet','Energy_mace']].mean(axis=1)
        # merged_df = merged_df.drop(columns=['Energy_1', 'Energy_2', 'Energy_3'])

        merged_df=merged_df.sort_values(by=["Mean_Energy"])

        merged_df.to_csv(path_results+"MajorityVote_"+formula+"_all.csv",index=False)
        # print("Majority vote is successfully determined.")

        def extract_symmetry(name):
            match = name.split("_")[2]
            return match

        merged_df["Sym"] = merged_df["CIF_Name"].apply(extract_symmetry)

        merged_df=merged_df.drop_duplicates(subset=["Sym"])
        merged_df = merged_df.drop(columns=["Sym"])
        merged_df.to_csv(path_results+"MajorityVote_"+formula+"_best.csv",index=False)
        
    if N_model==2:
        df_m3gnet=pd.read_csv(path_m3gnet+"M3Gnet_"+formula+"_all.csv")
        df_mace=pd.read_csv(path_mace+"MACE_"+formula+"_all.csv")

        # merged_df=pd.merge(df_alignn, df_m3gnet, on='CIF_Name', how='inner')
        merged_df = pd.merge(df_mace, df_m3gnet, on=['CIF_Name', 'Neigh_Order'], suffixes=('_mace', '_m3gnet'))

        merged_df['Mean_Energy'] = merged_df[['Energy_m3gnet','Energy_mace']].mean(axis=1)
        # merged_df = merged_df.drop(columns=['Energy_1', 'Energy_2', 'Energy_3'])

        merged_df=merged_df.sort_values(by=["Mean_Energy"])

        merged_df.to_csv(path_results+"MajorityVote_"+formula+"_all.csv",index=False)
        # print("Majority vote is successfully determined.")

        def extract_symmetry(name):
            match = name.split("_")[2]
            return match

        merged_df["Sym"] = merged_df["CIF_Name"].apply(extract_symmetry)

        merged_df=merged_df.drop_duplicates(subset=["Sym"])
        merged_df = merged_df.drop(columns=["Sym"])
        merged_df.to_csv(path_results+"MajorityVote_"+formula+"_best.csv",index=False)
    print("Majority vote is successfully determined.")

















# def m3gnet_energy_test():
#     # import warnings
#     from m3gnet.models import M3GNet
#     from pymatgen.core import Lattice, Structure

#     # suppress TensorFlow warnings
#     # for category in (UserWarning, DeprecationWarning):
#     #     warnings.filterwarnings("ignore", category=category, module="tensorflow")

#     # 1. Build your (possibly unrelaxed) structure:
#     #    Here we use Mo in a slightly stretched cell.
#     mo = Structure(
#         Lattice.cubic(3.3), 
#         ["Mo", "Mo"], 
#         [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
#     )

#     # 2. Load the pretrained M3GNet model:
#     #    This is a shortcut for M3GNet.from_dir(...) that downloads weights if needed.
#     model = M3GNet.load()  # defaults to "MP-2021.2.8-EFS" :contentReference[oaicite:1]{index=1}

#     # 3. Predict the (intensive) energy:
#     #    Returns a tf.Tensor of shape () giving eV/atom.
#     energy_per_atom = model.predict_structure(mo)

#     # 4. (Optional) convert to a Python float:
#     energy_per_atom = float(energy_per_atom)

#     print(f"Predicted energy (no relaxation) = {energy_per_atom:.4f} eV/atom")

# def megnet_calc_old(path,path_results):
#     import warnings

#     from m3gnet.models import Relaxer
#     from pymatgen.core import Lattice, Structure

#     import os
#     import pandas as pd

#     for category in (UserWarning, DeprecationWarning):
#         warnings.filterwarnings("ignore", category=category, module="tensorflow")

#     from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
#     # path="./Unknown_Systems_Data/Negative/"
    
#     print("\nM3GNET is activated...\n")
#     if not os.path.exists(path_results):
#         os.makedirs(path_results)
#     comp_list=os.listdir(path)
#     for comp in comp_list:
#         if not "output_" in comp:
#             continue
#         print(comp)
#         struct_list_init=[]
#         neigh_list=os.listdir(path+comp)
#         if(neigh_list==[]):
#             continue
#         for neigh in neigh_list:
#             if "Neigh" not in neigh:
#                 continue
#             # if neigh!="1_Neigh":
#             #     continue
#             sym_list=os.listdir(path+comp+"/"+neigh)
#             if(sym_list==[]):
#                 continue
#             for sym in sym_list:
#                 proto_list=os.listdir(path+comp+"/"+neigh+"/"+sym)
#                 if(proto_list==[]):
#                     continue
#                 for proto in proto_list:
#                     target_path=path+comp+"/"+neigh+"/"+sym+"/"+proto
#                     # print(comp,neigh,sym,proto)
#                     struct=Structure.from_file(target_path)
#                     struct_list_init.append([struct,comp,neigh,sym,proto])
#                     # print(proto)
#         print("\nTotal number of systems:",len(struct_list_init))
#         print("------------------------")
#         print("Relaxer starts...")
        
#         final_structure_list=[]
#         relaxer=Relaxer()
#         for struct in struct_list_init:
#             print(struct[-1])
#             relax_results = relaxer.relax(struct[0], verbose=False)

#             final_structure=relax_results['final_structure']
#             final_energy_per_atom=float(relax_results['trajectory'].energies[-1] / len(struct[0]))
#             final_structure_list.append({"CIF_Name":struct[4].replace(".cif",""),"Sym":struct[3],"Neigh_Order":struct[2],"Energy":final_energy_per_atom,"Target":struct[1].replace("output_",""),"Struc_Init":struct[0],"Struc_Final":final_structure,})
#             if not os.path.exists(path_results+struct[4]+"_opt.cif"):
#                 final_structure.to_file(path_results+struct[4]+"_opt.cif")
#             # print(f"Relaxed lattice parameter is {final_structure.lattice.abc[0]:.3f} Å")
#             # print(f"Final energy is {final_energy_per_atom:.3f} eV/atom")
#         df=pd.DataFrame(final_structure_list)
#         df['Energy'] = df['Energy'].round(3)
#         df=df.sort_values(by=["Energy"])
#         df=df.drop_duplicates(subset=["Sym"])
#         df[["CIF_Name","Sym","Neigh_Order","Energy","Target"]].to_csv(path_results+final_structure_list[0]["Target"].replace("output_","")+'.csv', index=False)
#         print("Done")
#         print("------------------------")
#     print("Terminated Sucessfully!")