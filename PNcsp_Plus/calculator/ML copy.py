def ALIGNN_calc(formula,path,path_results,relax):
    from alignn.ff.ff import AlignnAtomwiseCalculator, default_path
    from jarvis.core.atoms import Atoms, ase_to_atoms
    from ase.constraints import ExpCellFilter
    from ase.optimize.fire import FIRE
    import matplotlib.pyplot as plt
    import numpy as np

    import pandas as pd
    import os
    from pymatgen.core import Lattice, Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from ase.io import write



    print("\nALIGNN is activated...\n")
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
                struct = Atoms.from_cif(target_path)
                struct_list_init.append([struct,formula,neigh,sym,proto])
                # print(proto)
    print("\nTotal number of systems:",len(struct_list_init))
    print("------------------------")
    print("Relaxer starts:")
    
    final_structure_list=[]

    # Relaxer function
    def general_relaxer(atoms, calculator, fmax=0.05, steps=150, relax=relax, CellFiter=True):
        print(f"Number of atoms in CIF: {atoms.num_atoms}")
        ase_atoms = atoms.ase_converter()
        print("ATOMS")
        print(atoms)
        print("ASE_ATOMS")
        print(ase_atoms.positions)
        ase_atoms.calc = calculator
        if not relax:
            return ase_atoms, ase_atoms.get_potential_energy()
        
        if(CellFiter==True):
            ase_atoms = ExpCellFilter(ase_atoms)
            dyn = FIRE(ase_atoms)
            dyn.run(fmax=fmax, steps=steps)
            # return ase_to_atoms(ase_atoms.atoms), ase_atoms.get_potential_energy()
            return ase_atoms.atoms,  ase_atoms.atoms.get_potential_energy()
        else:
            dyn = FIRE(ase_atoms)
            dyn.run(fmax=fmax, steps=steps)
            # return ase_to_atoms(ase_atoms.atoms), ase_atoms.get_potential_energy()
            return ase_atoms, ase_atoms.get_potential_energy()

    
    # ALIGNN Calculator
    model_path = default_path()
    # model_path = "/home/cem/alignn/alignn/ff/jv_formation_energy_peratom_alignn"
    calc = AlignnAtomwiseCalculator(path=model_path)

    # # Load CIF
    # atoms = Atoms.from_cif("./Si2.cif")

    # Main logic
    for struct in struct_list_init:
        # print(struct[-1])
        # print("Initial structure:\n", struct[0])
        # atoms = general_relaxer(atoms=struct, calculator=calc)
        atom_num=struct[0].num_atoms
        final_structure, final_energy = general_relaxer(atoms=struct[0], calculator=calc,CellFiter=True, relax=relax)
        # print("Initial structure:\n", struct[0])
        # print("Relaxed structure:\n", final_structure)
        print("Energy: ",final_energy/atom_num," eV/atom")
        final_structure_list.append({"Neigh_Name":struct[4].replace(".cif",""),"Sym":struct[3],"Neigh_Order":struct[2],"Energy":final_energy/atom_num,"Target":struct[1],"Struc_Init":struct[0],"Struc_Final":final_structure})
        print(type(final_structure))

        if not os.path.exists(path_results+struct[4]+"_opt.cif"):
            # final_structure.to_file(path_results+struct[4]+"_opt.cif")
            write(path_results+struct[4]+"_opt.cif", final_structure)
    
    df=pd.DataFrame(final_structure_list)
    # df['Energy'] = pd.to_numeric(df['Energy'], errors='coerce').round(3)
    # df['Energy'] = df['Energy'].round(3)
    df=df.sort_values(by=["Energy"])
    df[["Neigh_Name","Sym","Neigh_Order","Energy","Target"]].to_csv(path_results+final_structure_list[0]["Target"]+'_all.csv', index=False)

    df=df.drop_duplicates(subset=["Sym"])
    df[["Neigh_Name","Sym","Neigh_Order","Energy","Target"]].to_csv(path_results+final_structure_list[0]["Target"]+'_best.csv', index=False)
    print("Done")
    print("------------------------")
    print("Terminated Sucessfully!")


def megnet_calc(formula,path,path_results):
    import warnings

    from m3gnet.models import Relaxer
    from pymatgen.core import Lattice, Structure

    import os
    import pandas as pd

    for category in (UserWarning, DeprecationWarning):
        warnings.filterwarnings("ignore", category=category, module="tensorflow")

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
                struct=Structure.from_file(target_path)
                struct_list_init.append([struct,formula,neigh,sym,proto])
                # print(proto)
    print("\nTotal number of systems:",len(struct_list_init))
    print("------------------------")
    print("Relaxer starts:")
    
    final_structure_list=[]
    relaxer=Relaxer()
    for struct in struct_list_init:
        print(struct[-1])
        relax_results = relaxer.relax(struct[0], verbose=False)

        final_structure=relax_results['final_structure']
        final_energy_per_atom=float(relax_results['trajectory'].energies[-1] / len(struct[0]))
        final_structure_list.append({"Neigh_Name":struct[4].replace(".cif",""),"Sym":struct[3],"Neigh_Order":struct[2],"Energy":final_energy_per_atom,"Target":struct[1],"Struc_Init":struct[0],"Struc_Final":final_structure,})
        if not os.path.exists(path_results+struct[4]+"_opt.cif"):
            final_structure.to_file(path_results+struct[4]+"_opt.cif")
        # print(f"Relaxed lattice parameter is {final_structure.lattice.abc[0]:.3f} Å")
        # print(f"Final energy is {final_energy_per_atom:.3f} eV/atom")
    df=pd.DataFrame(final_structure_list)
    df['Energy'] = df['Energy'].round(3)
    df=df.sort_values(by=["Energy"])
    df[["Neigh_Name","Sym","Neigh_Order","Energy","Target"]].to_csv(path_results+final_structure_list[0]["Target"]+'_all.csv', index=False)

    df=df.drop_duplicates(subset=["Sym"])
    df[["Neigh_Name","Sym","Neigh_Order","Energy","Target"]].to_csv(path_results+final_structure_list[0]["Target"]+'_best.csv', index=False)
    print("Done")
    print("------------------------")
    print("Terminated Sucessfully!")


def megnet_calc_old(path,path_results):
    import warnings

    from m3gnet.models import Relaxer
    from pymatgen.core import Lattice, Structure

    import os
    import pandas as pd

    for category in (UserWarning, DeprecationWarning):
        warnings.filterwarnings("ignore", category=category, module="tensorflow")

    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    # path="./Unknown_Systems_Data/Negative/"
    
    print("\nM3GNET is activated...\n")
    if not os.path.exists(path_results):
        os.makedirs(path_results)
    comp_list=os.listdir(path)
    for comp in comp_list:
        if not "output_" in comp:
            continue
        print(comp)
        struct_list_init=[]
        neigh_list=os.listdir(path+comp)
        if(neigh_list==[]):
            continue
        for neigh in neigh_list:
            if "Neigh" not in neigh:
                continue
            # if neigh!="1_Neigh":
            #     continue
            sym_list=os.listdir(path+comp+"/"+neigh)
            if(sym_list==[]):
                continue
            for sym in sym_list:
                proto_list=os.listdir(path+comp+"/"+neigh+"/"+sym)
                if(proto_list==[]):
                    continue
                for proto in proto_list:
                    target_path=path+comp+"/"+neigh+"/"+sym+"/"+proto
                    # print(comp,neigh,sym,proto)
                    struct=Structure.from_file(target_path)
                    struct_list_init.append([struct,comp,neigh,sym,proto])
                    # print(proto)
        print("\nTotal number of systems:",len(struct_list_init))
        print("------------------------")
        print("Relaxer starts:")
        
        final_structure_list=[]
        relaxer=Relaxer()
        for struct in struct_list_init:
            print(struct[-1])
            relax_results = relaxer.relax(struct[0], verbose=False)

            final_structure=relax_results['final_structure']
            final_energy_per_atom=float(relax_results['trajectory'].energies[-1] / len(struct[0]))
            final_structure_list.append({"Neigh_Name":struct[4].replace(".cif",""),"Sym":struct[3],"Neigh_Order":struct[2],"Energy":final_energy_per_atom,"Target":struct[1].replace("output_",""),"Struc_Init":struct[0],"Struc_Final":final_structure,})
            if not os.path.exists(path_results+struct[4]+"_opt.cif"):
                final_structure.to_file(path_results+struct[4]+"_opt.cif")
            # print(f"Relaxed lattice parameter is {final_structure.lattice.abc[0]:.3f} Å")
            # print(f"Final energy is {final_energy_per_atom:.3f} eV/atom")
        df=pd.DataFrame(final_structure_list)
        df['Energy'] = df['Energy'].round(3)
        df=df.sort_values(by=["Energy"])
        df=df.drop_duplicates(subset=["Sym"])
        df[["Neigh_Name","Sym","Neigh_Order","Energy","Target"]].to_csv(path_results+final_structure_list[0]["Target"].replace("output_","")+'.csv', index=False)
        print("Done")
        print("------------------------")
    print("Terminated Sucessfully!")