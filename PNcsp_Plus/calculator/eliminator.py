import os
import pandas as pd
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
import numpy as np

def eliminator(formula,path):
    struct_list_init=[]
    neigh_list=os.listdir(path)
    if(neigh_list==[]):
        exit(0)
    for neigh in neigh_list:
        if "Neigh" not in neigh:
            continue

        sym_list=os.listdir(path+neigh)
        if(sym_list==[]):
            continue
        for sym in sym_list:
            proto_list=os.listdir(path+neigh+"/"+sym)
            if(proto_list==[]):
                continue
            for proto in proto_list:
                target_path=path+neigh+"/"+sym+"/"+proto

                struct=Structure.from_file(target_path)

                struct_list_init.append([formula,neigh,sym,proto,target_path,struct])
                
    return struct_list_init,sym_list

path0="/home/cem/PNcsp_Plus/PNcsp/DATA/DATA_180_4NN_eliminated/"

formula_list=[x.replace("output_","") for x in os.listdir(path0) if "output_" in x]

for formula in formula_list:  
    print(formula)  
    path=path0+"output_"+formula+"/"
    struct_list_init,sym_list=eliminator(formula,path)

    df=pd.DataFrame(struct_list_init, columns=["formula","neigh","sym","proto","target_path","struct"])

    sm = StructureMatcher()
    for sym in sym_list:
        df_list=list(df[df["sym"]==sym].sort_values(by="neigh",ascending=False).values)
        for i in range(len(df_list) - 1, -1, -1):
            # if not os.path.exists(df_list[i][4]):
            #     continue
            for j in range(i - 1, -1, -1):
                if not os.path.exists(df_list[j][4]):
                    continue
                my_struct1=df_list[i][5]
                my_struct2=df_list[j][5]
                fit_result=sm.fit(my_struct1, my_struct2)
                if(fit_result==True):
                    os.remove(df_list[j][4])
                    # del df_list[j]
                    # break