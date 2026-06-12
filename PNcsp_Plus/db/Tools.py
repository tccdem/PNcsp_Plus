import os 
import pandas as pd
def structure_collector(path):
    from pymatgen.core import Structure
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

    struc_df=pd.DataFrame(struct_list_init,columns=["CIF_Name","Neigh_Order", "sym", "struct"])

    return struc_df

def table_matcher(path,struc_df,GNN):

    dest_path=os.path.join(path,"Calc_report")
    if not os.path.exists(dest_path):
        pass
    else:
        folder_list = [name for name in os.listdir(dest_path) if os.path.isdir(os.path.join(dest_path, name))]
        for folder in folder_list:
            if GNN not in folder:
                continue
            csv_path=os.path.join(dest_path,folder)
            csv_list=os.listdir(csv_path)
            for csv_file in csv_list:
                if "_all" in csv_file:
                    GNN_df=pd.read_csv(os.path.join(csv_path,csv_file))
                
            common_cols = [col for col in struc_df.columns if col in GNN_df.columns]
            extra_cols = [col for col in GNN_df.columns if col not in struc_df.columns]

            struc_df_matched = struc_df.merge(
                GNN_df[common_cols + extra_cols].drop_duplicates(subset=common_cols),
                on=common_cols,
                how="left"
            )

            struc_df_matched=struc_df_matched.sort_values(by="Energy")
    
    return struc_df_matched

def data_reduction(path,GNN):
    from pymatgen.analysis.structure_matcher import StructureMatcher
    struc_df=structure_collector(path)
    struc_df_matched=table_matcher(path,struc_df,GNN)

    sm = StructureMatcher()
    all_rows = []
    for sym, group in struc_df_matched.groupby("sym"):
        # group = group.sort_values(by=["Neigh_Order", "Energy"],ascending=[False, True]).reset_index(drop=True)
        representatives = []
        duplicate_of = []
        for i, row in group.iterrows():
            current_struct = row["struct"]

            matched = False
            matched_rep_index = None

            for rep_idx, rep_struct in representatives:
                if sm.fit(current_struct, rep_struct):
                    matched = True
                    matched_rep_index = rep_idx
                    # matched_rep_index = row["CIF_Name"]
                    break

            if matched:
                duplicate_of.append(matched_rep_index)
            else:
                representatives.append((row["CIF_Name"], current_struct))
                duplicate_of.append(None)
        group["duplicate_of"] = duplicate_of
        group["is_unique"] = group["duplicate_of"].isna()
    
        all_rows.append(group)

    my_df=pd.concat(all_rows, ignore_index=True)
    my_df.drop(["struct"],axis=1,inplace=True)
    my_df = my_df[my_df["duplicate_of"].isna()]
    my_df.drop(["duplicate_of","is_unique"],axis=1,inplace=True)
    my_df = my_df.sort_values(by="Energy")
    return my_df
    