import os
import re
import pandas as pd
def offline_order(Comp_list):

    Comp_list_ordered=[]
    for i in range(len(Comp_list)):
        comp_id=[x for x in re.split(r'(?=[A-Z])', Comp_list[i]) if x]
        comp_id=sorted(comp_id)
        comp_query=""
        for j in range(len(comp_id)):
            # if not re.search(r'\d', comp_id[i]):
            #     comp_id[i]=comp_id[i]+"1"
            comp_query+=comp_id[j]
            if (j!=len(comp_id)-1):
                comp_query+=" "
        Comp_list_ordered.append(comp_query)
    print(Comp_list_ordered)
    return Comp_list_ordered

def get_OQMD_data(formula,path):

    from ase.spacegroup import crystal
    from ase.io import write
    from itertools import product
    from math import gcd
    from functools import reduce
    import shutil
    import qmpy

    def write_CIF(elem_list, dest_path, formula,spacegroup_num, ind_ext):
        struct = crystal(elem_list, site_list, cell=cell, size=(1, 1, 1))
        filename = f"{dest_path}{formula}_sym{spacegroup_num}_OQMD_{ind_ext+1}.cif"
        write(filename, struct)

    def reduce_list(nums):
        g = reduce(gcd, nums)
        return [n // g for n in nums]

    
    def parse_formula_counts(formula):
        matches = re.findall(r'[A-Z][a-z]*|\d+', formula)
        return {matches[i]: int(matches[i+1]) for i in range(0, len(matches), 2)}

    All_list=[]
    Comp_ordered=offline_order([formula])[0]
    raw_data = qmpy.Entry.objects.filter(composition_id=Comp_ordered)
    if(len(raw_data)==0):
        print(Comp_ordered,"--> no structure")
        print("Target compound can not be found in OQMD!!")
        return None

    target_data={"data":[]}
    for k in range(len(raw_data)):
        comp=raw_data[k]
        delta_e=comp.energy
        target_data["data"].append({"name":comp.name, "spacegroup":comp.spacegroup, "cell":comp.structure.cell, "sites":comp.structure.atoms,"delta_e":delta_e})

    dest_path = f"{path}Known_Structures/"

    if not os.path.exists(dest_path):
            os.makedirs(dest_path)

    for ind_ext, entry in enumerate(target_data['data']):
        name = entry['name']
        # spacegroup_sym = entry['spacegroup'].symbol
        spacegroup_num = entry['spacegroup'].number
        cell = entry['cell']
        sites = entry['sites']

        # Collect element list and site list
        site_list, elem_list, elem_list_org = [], [], []
        for s in sites:
            elem, coords = str(s).split(' @ ')
            site_list.append(tuple(coords.split()))
            elem_list.append(elem)      

        write_CIF(elem_list, dest_path, formula,spacegroup_num, ind_ext)

    print("OQMD data is collected!")
    return target_data

def get_MP_data(formula,path):
    from mp_api.client import MPRester

    def write_CIF(struct, dest_path, formula,spacegroup_num, ind_ext):
        filename = f"{dest_path}{formula}_sym{spacegroup_num}_MP_{ind_ext+1}.cif"
        # write(filename, struct)
        struct.to(fmt="cif", filename=filename)

    Key_path="./dev/MP_API_KEY"
    with open(Key_path, "r") as f:
        MP_API_KEY = f.readline().strip()

    with MPRester(MP_API_KEY) as mpr:
        fields = [
        "formula_pretty",
        "composition_reduced",
        # "composition",
        # "elements",
        "symmetry",
        "structure",
        "energy_per_atom",
        "formation_energy_per_atom"
        ]
        target_data=mpr.summary.search(formula=formula,fields=fields)

    dest_path = f"{path}Known_Structures/"

    if not os.path.exists(dest_path):
            os.makedirs(dest_path)

    for ind_ext, entry in enumerate(target_data):
        name = entry.formula_pretty
        structure = entry.structure
        # spacegroup_sym = entry.symmetry.symbol
        spacegroup_num = entry.symmetry.number

        write_CIF(structure, dest_path, formula,spacegroup_num, ind_ext)
    
    print("MP data is collected!")
    return target_data

def compare_structures(formula,path):

    from pymatgen.analysis.structure_matcher import StructureMatcher
    from pymatgen.core import Structure

    def ref_collector(path):
        structure_ref = []
        known_path=os.path.join(path, "Known_Structures")
        struct_list=os.listdir(known_path)
        for struct_name in struct_list:
            target_path=os.path.join(known_path,struct_name)
            sym_num=struct_name.split("sym")[1].split("_")[0]
            try:
                struct = Structure.from_file(target_path)
                structure_ref.append([struct_name.replace(".cif",""),sym_num,struct])
            except Exception as e:
                print(f"Could not read {target_path}: {e}")
        return structure_ref
    

    def predicted_collector(path):
        struct_list_init = []

        neigh_list = [x for x in os.listdir(path) if "Neigh" in x]
        if not neigh_list:
            return struct_list_init

        for neigh in neigh_list:
            neigh_path = os.path.join(path, neigh)
            if not os.path.isdir(neigh_path):
                continue

            sym_list = os.listdir(neigh_path)
            if not sym_list:
                continue

            for sym in sym_list:
                sym_path = os.path.join(neigh_path, sym)
                if not os.path.isdir(sym_path):
                    continue

                proto_list = os.listdir(sym_path)
                if not proto_list:
                    continue

                for struct_name in proto_list:
                    target_path = os.path.join(sym_path, struct_name)

                    if not os.path.isfile(target_path):
                        continue


                    try:
                        struct = Structure.from_file(target_path)
                        struct_list_init.append([struct_name.replace(".cif",""),neigh,int(sym.replace("sym","")), struct])
                    except Exception as e:
                        print(f"Could not read {target_path}: {e}")

        return struct_list_init
    
    ref_df=pd.DataFrame(ref_collector(path),columns=["CIF_Name","sym", "struct"])
    my_df=pd.DataFrame(predicted_collector(path),columns=["CIF_Name","Neigh_Order", "sym", "struct"])
    sm = StructureMatcher()

    unique_list_sym=[]
    unique_list_struc=[]

    ref_df["sym"]
    sym_ref_list=list(set(ref_df["sym"].astype(int).values))
    for i, my_row in my_df.iterrows():
        my_struct = my_row["struct"]
        is_struc_unique=True
        is_sym_unique=True
        if(my_row["sym"] in sym_ref_list):
            is_sym_unique=False

        for j, ref_row in ref_df.iterrows():
            ref_struct = ref_row["struct"]
            if sm.fit(my_struct, ref_struct):
                is_struc_unique=False
                break
        unique_list_sym.append(is_sym_unique)
        unique_list_struc.append(is_struc_unique)
    my_df["is_new_sym"]=unique_list_sym
    my_df["is_new_struc"]=unique_list_struc

    print("New structures are successfully detected!!")
    my_df.drop(["struct"],axis=1,inplace=True)
    dest_path=os.path.join(path,"Calc_report")
    my_df_new=my_df[(my_df["is_new_struc"]==True)].copy()
    my_df_new.drop(["is_new_struc"],axis=1,inplace=True)
    return my_df,my_df_new

def matcher(formula,path,struc_df,new_struc_df):
    dest_path=os.path.join(path,"Calc_report")
    if not os.path.exists(dest_path):
        os.makedirs(dest_path)
        struc_df.to_csv(os.path.join(dest_path,"Similarity_Report.csv"),index=False)
        new_struc_df.to_csv(os.path.join(dest_path,"New_Structures.csv"),index=False)
    else:
        struc_df.to_csv(os.path.join(dest_path,"Similarity_Report.csv"),index=False)
        folder_list = [name for name in os.listdir(dest_path) if os.path.isdir(os.path.join(dest_path, name))]
        for folder in folder_list:
            csv_path=os.path.join(dest_path,folder)
            csv_list=os.listdir(csv_path)
            for csv_file in csv_list:
                if "_all" in csv_file:
                    GNN_df=pd.read_csv(os.path.join(csv_path,csv_file))
                
            common_cols = [col for col in new_struc_df.columns if col in GNN_df.columns]
            extra_cols = [col for col in GNN_df.columns if col not in new_struc_df.columns]

            new_struc_df = new_struc_df.merge(
                GNN_df[common_cols + extra_cols].drop_duplicates(subset=common_cols),
                on=common_cols,
                how="left"
            )

            struc_df = struc_df.merge(
                GNN_df[common_cols + extra_cols].drop_duplicates(subset=common_cols),
                on=common_cols,
                how="left"
            )

            new_struc_df=new_struc_df.sort_values(by="Energy")
            new_struc_df.to_csv(os.path.join(dest_path,folder,folder+"_"+formula+"_newstruc_all.csv"),index=False)

            new_struc_df=new_struc_df.drop_duplicates(subset=["sym"])
            new_struc_df.to_csv(os.path.join(dest_path,folder,folder+"_"+formula+"_newstruc_best.csv"),index=False)

            struc_df=struc_df.sort_values(by="Energy")
            struc_df.to_csv(os.path.join(dest_path,folder,folder+"_"+formula+"_all.csv"),index=False)

            struc_df=struc_df.drop_duplicates(subset=["sym"])
            struc_df.to_csv(os.path.join(dest_path,folder,folder+"_"+formula+"_best.csv"),index=False)


def find_unique_data(formula,path):

    get_OQMD_data(formula,path)
    get_MP_data(formula,path)
    struc_df,new_struc_df=compare_structures(formula,path)
    matcher(formula,path,struc_df,new_struc_df)

    