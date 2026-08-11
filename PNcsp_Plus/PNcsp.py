#!/usr/bin/env python

###############################################################
#                                                             #
#                          PNcsp+                             #
#                                                             #
###############################################################

import numpy as np
import re
import os
import itertools
import time
import sys
import argparse
from pathlib import Path
from pymatgen.core import Composition

from PNcsp_Plus.calculator import ML
from PNcsp_Plus.db import DBsearch

# Get the directory where PNcsp.py actually lives
PACKAGE_DIR = Path(__file__).resolve().parent

def blockPrint():
    sys.stdout = open(os.devnull, 'w')
def enablePrint():
    sys.stdout = sys.__stdout__

def get_Symbol(PN):
    csv_path = PACKAGE_DIR / "db" / "data" / "Z_PN_Elem_extended_MD.csv"
    with open(csv_path,'r') as file:
        Comp_list={}
        for line in file:
            if 'PN' in line:
                continue
            info=line.strip().replace('\ufeff','').split(',')
            Comp_list[info[0]]=info[1]
    return(Comp_list[str(PN)])

def get_PN(Symb):

    csv_path = PACKAGE_DIR / "db" / "data" / "Z_PN_Elem_extended_MD.csv"
    with open(csv_path,'r') as file:
        Comp_list={}
        for line in file:
            if 'PN' in line:
                continue
            info=line.strip().replace('\ufeff','').split(',')
            Comp_list[info[1]]=int(info[0])
    # print(Comp_list)
    return(Comp_list[Symb])

def separate_elems(formula):
  """Separates a chemical formula using regular expressions"""
  matches = re.findall(r'[A-Z][a-z]*|\d+', formula)
  elems = [match for match in matches if match.isalpha()]
  counts = [int(match) for match in matches if match.isdigit()]
  return elems, counts

def convert_formula(formula):
    elems, counts=separate_elems(formula)
    PNs=[]
    for Symb in elems:
        PNs.append(get_PN(Symb))
    print("Separated Formula:",elems)
    print("PNs:",PNs)
    print("Ratios:",counts)
    return elems, counts, PNs

def dist_classifier(PNs,res):
    dists=[]
    for i in range(len(res)):
        dist_local=[]
        for j in range(len(PNs)):
            # for k in range(len(res[i])):
            dist_local.append(abs(PNs[j]-res[i][j]))
        dist_max=max(dist_local)
        dists.append(dist_max)
    return dists

def get_Neig(formula, N_neig):
    elems, counts, PNs=convert_formula(formula)

    # Generate neigh PN list (exclude the phase map borders)
    PN_new_all=[]
    for i in range(len(PNs)):
        PN_new=[]
        for j in np.arange(-1*N_neig,N_neig+1):
            if((PNs[i]+j>0) and (PNs[i]+j<119)):
                PN_new.append(PNs[i]+j)

        PN_new_all.append(PN_new)
    print("PN_new_all: ",PN_new_all)

    # In case constituent elemens are too close 
    for i in range(len(PNs)-1):
        for j in range(i+1,len(PNs)):
            if(abs(PNs[i]-PNs[j])<=2*N_neig):
                print("\nWARNING: The neighbors of the ",elems[i],"and",elems[j] ,"overlap. This is not a problem, but it is advisable to exercise caution when examining the output CIFs.")
    # New Exchange Look-up Table
    exchange_dict={}
    for i in range(len(PN_new_all)):
        for j in range(len(PN_new_all[i])):
            key=get_Symbol(PN_new_all[i][j])
            if exchange_dict.get(key) is not None:
                exchange_dict[key].append(elems[i])
            else:
                exchange_dict[key]=[elems[i]]
    print("exchange_dict: ",exchange_dict)

    # Exchange Look-up Table for PN
    exchange_dict_PN={}
    for i in range(len(PN_new_all)):
        for j in range(len(PN_new_all[i])):
            key=PN_new_all[i][j]
            if exchange_dict_PN.get(key) is not None:
                exchange_dict_PN[key].append(PNs[i])
            else:
                exchange_dict_PN[key]=[PNs[i]]

    res=np.array(list(itertools.product(*PN_new_all)))

    drop_list=[]
    for i in range(len(res)):
        # Drop original formula
        if(np.array_equal(sorted(res[i]), sorted(PNs))):
            drop_list.append(i)
        # Drop elemental dublication
        if(len(set(res[i]))<len(res[i])):
            drop_list.append(i)

    res=np.delete(res, drop_list,axis=0)

    # Distance Detector
    dist_list=dist_classifier(PNs,res)

    # Order according to neigh distance
    res_ordered=[]
    dist_list_ordered=[]
    for dst in range(1,N_neig+1):
        for i in range(len(dist_list)):
            if (dist_list[i]==dst):
                res_ordered.append(list(res[i]))
                dist_list_ordered.append(dist_list[i])

    Symbol_list=[]
    for i in range(len(res_ordered)):
        Symbols=[]
        for j in range(len(res_ordered[i])):
            Symbols.append(get_Symbol(res_ordered[i][j]))
        Neig_formula=""    
        for k in range(len(Symbols)):
            Neig_formula+=Symbols[k]+str(counts[k])
        Symbol_list.append(Neig_formula)

    for i in range(len(dist_list_ordered)):
        if(dist_list_ordered[i]==N_neig):
            res_ordered=res_ordered[i:]
            dist_list_ordered=dist_list_ordered[i:]
            Symbol_list=Symbol_list[i:]
            break

    return Symbol_list,dist_list_ordered,exchange_dict

def categorize(N_neig,formula,data_path):
    import shutil  
    import os
    
    path=data_path+"/output_"+formula+"/"+str(N_neig)+"_Neigh/"
    cifs=[f for f in os.listdir(path) if ".cif" in f]
    print("Total Number of CIF files:",len(cifs))
    for cif in cifs:
        source_path=path+cif

        # sym=cif.split("_")[2]
        sym="sym"+cif.split("sym")[1].split("_")[0]
        dest_path=path+sym+"/"

        if not os.path.exists(dest_path):
            os.makedirs(dest_path)

        dest = shutil.move(source_path, dest_path)  

def show_config(formula,N_neig,E_filter,timer,online,calculator,database,BlockSearch,Relaxer,data_path,CheckNew,top_n,top_c,Reverse,StoreData):
    print("\nProgram Configuration")
    print("---------------------")
    print("Query formula:\t",formula,"\nNeighbor order:\t",N_neig,
          "\nEnergy filter:\t",E_filter,"\nSleep timer:\t",timer,
          "\nOnline: \t",online,"\nCalculator: \t",calculator,"\nData source: \t",database,
          "\nBlockSearch: \t",BlockSearch,"\nRelaxer: \t",Relaxer,
          "\nOutputDir: \t",data_path,"\nCheckNew: \t",CheckNew,"\nTop n:  \t",top_n,"\nTop c:  \t",top_c,"\nReverse: \t",Reverse,"\nStore Data: \t",StoreData)
    print("---------------------\n")

def order_formula(comp):
    import re
    
    parts = re.findall(r"([A-Z][a-z]?)(\d*)", comp.strip())
    comp_w1=""
    for el, num in parts:
        comp_w1+=el
        comp_w1+=num if num else "1"

    comp_decomp=[x for x in re.split(r'(?=[A-Z])', comp_w1) if x]
    # comp_decomp=sorted(comp_decomp)
    # comp_ordered=""
    # for j in range(len(comp_decomp)):
    #     comp_ordered+=comp_decomp[j]

    # if len(re.sub(r"\d+", "", comp_ordered))!=len(re.sub(r"\d+", "", comp)):
    #     print(f"Error: Invalid formula: '{comp}'. Elements start with capital letter and followed by ration. For ex: Na1Cl1, Al1Fe2")
    #     return ""
    return comp_w1  

def run(formula,N_neig=1,E_filter=0,time_sleep=None,online="False",calculator=None,data_path=".",
        BlockSearch=False,Relax=False,database="OQMD",CheckNew=False,top_n=None,
        top_c=None,project=None,dim=None,ReduceData=False,Reverse=False,StoreData=False):
    
    formula=order_formula(str(Composition(formula)).replace(" ", ""))
    # formula=order_formula(formula)

    if not formula:
        return(None)

    show_config(formula=formula,N_neig=N_neig,E_filter=E_filter,timer=time_sleep,online=online,
                calculator=calculator,database=database,BlockSearch=BlockSearch,Relaxer=Relax,
                data_path=data_path,CheckNew=CheckNew,top_n=top_n,top_c=top_c,Reverse=Reverse,StoreData=StoreData)
    path=data_path+"/output_"+formula+"/"

    if(Reverse==True):
        res,neigh_list,exchange_dict=get_Neig(formula=formula,N_neig=N_neig)
        neigh_path=os.path.join(data_path,formula+"_neighborhood")
        os.makedirs(neigh_path,exist_ok=True)
        
        with open(os.path.join(neigh_path,"neighborhood_"+str(N_neig)+".txt"),"w") as f:
            f.write("System\tKnown_Structures\n")
            for ne in res:
                DBsearch.get_OQMD_data(ne,neigh_path+"/"+ne+"_")
                DBsearch.get_MP_data(ne,neigh_path+"/"+ne+"_")
                f.write(ne+"\t")
                pos_count=0
                neg_count=0
                if os.path.exists(neigh_path+"/"+ne+"_Known_Structures"):
                    for fname in os.listdir(neigh_path+"/"+ne+"_Known_Structures"):
                        if fname.endswith("pos.cif"):
                            pos_count += 1
                        elif fname.endswith("neg.cif"):
                            neg_count += 1
                f.write(" pos:"+str(pos_count)+" neg:"+str(neg_count)+"\n")
        return(None)

        
    if(BlockSearch!=True):
        for neig in range(1,N_neig+1):
            print("\n*** Neighbor:",neig,"***")
            res,neigh_list,exchange_dict=get_Neig(formula=formula,N_neig=neig)

            if(database=="OQMD"):
                if(online==True):
                    from PNcsp_Plus.db import OQMDonline
                    All_list=OQMDonline.get_data_OQMD(res,neigh_list,Energy_filter=E_filter,timer=time_sleep)
                    if not All_list:
                        continue
                        return(None)
                    OQMDonline.create_prototype_OQMD(All_list,exchange_dict,formula=formula,data_path=data_path)
                else:
                    from PNcsp_Plus.db import OQMDoffline
                    All_list=OQMDoffline.get_data_OQMD(res,neigh_list,Energy_filter=E_filter)
                    if not All_list:
                        continue
                        return(None)
                    OQMDoffline.create_prototype_OQMD(All_list,exchange_dict,formula=formula,data_path=data_path)
            elif(database=="MP"):
                from PNcsp_Plus.db import MPonline
                All_list=MPonline.get_data_MP(res,Energy_filter=E_filter)
                if not All_list:
                    continue
                    return(None)
                MPonline.create_prototype_MP(All_list, exchange_dict, formula=formula, neigh=N_neig,data_path=data_path)
            elif(database=="MPDS"):
                print("WARNING: MPDS implementation is under construction. Choose another data source and run again.")
                return(None)
            categorize(N_neig=neig,formula=formula,data_path=data_path)
            
        print("STRUCTURE SEARCH TERMINATED SUCCESFULLY!\n")
        
    
    if(calculator=="M3GNet"):
        # ML.M3GNet_calc("./","./output_"+formula+"/Calc_report/")
        path_results=data_path+"/output_"+formula+"/Calc_report/M3GNet/"
        if(Relax==True):
            ML_res=ML.M3GNet_calc(formula,path,path_results,relax=True)
            if not ML_res:
                return(None)
        else:
            ML_res=ML.M3GNet_calc(formula,path,path_results,relax=False)
            if not ML_res:
                return(None)

    elif(calculator=="ALIGNN"):
        path_results=data_path+"/output_"+formula+"/Calc_report/ALIGNN/"
        if(Relax==True):
            ML_res=ML.ALIGNN_calc(formula,path,path_results,relax=True)
            if not ML_res:
                return(None)
        else:
            ML_res=ML.ALIGNN_calc(formula,path,path_results,relax=False)
            if not ML_res:
                return(None)
            
    elif(calculator=="MACE"):
        path_results=data_path+"/output_"+formula+"/Calc_report/MACE/"
        if(Relax==True):
            ML_res=ML.MACE_calc(formula,path,path_results,relax=True)
            if not ML_res:
                return(None)
        else:
            ML_res=ML.MACE_calc(formula,path,path_results,relax=False)
            if not ML_res:
                return(None)
    elif(calculator=="ensemble"):
        path_results=data_path+"/output_"+formula+"/Calc_report/ensemble/"
        if(Relax==True):
            print("\nERROR: You can not run Majority Vote with relaxation.\n")
        else:
            path_alignn=data_path+"/output_"+formula+"/Calc_report/ALIGNN/"
            ML_res=ML.ALIGNN_calc(formula,path,path_alignn,relax=False)
            if not ML_res:
                return(None)
            
            path_mace=data_path+"/output_"+formula+"/Calc_report/MACE/"
            ML.MACE_calc(formula,path,path_mace,relax=False)

            path_m3gnet=data_path+"/output_"+formula+"/Calc_report/M3GNet/"
            ML.M3GNet_calc(formula,path,path_m3gnet,relax=False)
            
            N_model=2 # 2 Model: m3gnet, mace 3 Model: alignn, m3gnet, mace
            ML.ensemble_vote(formula,path_alignn,path_mace,path_m3gnet,path_results,N_model)
    else:
        pass

    if(ReduceData==True): 
        from PNcsp_Plus.db import Tools
        import pandas as pd

        GNN_path=data_path+"/output_"+formula+"/Calc_report/"

        if not os.path.exists(GNN_path):
            print("Warning: GNN evaluation is missing. Run again with the flag -calc <model_name>!")
            return(None)
        Model_name_list=os.listdir(GNN_path)

        if "ensemble" in Model_name_list:
            my_df=Tools.data_reduction(path,"ensemble")
            dest_path_csv=os.path.join(path,"Calc_report/ensemble")
            dest_path_all=os.path.join(dest_path_csv, "ensemble_"+formula+"_all.csv")
            dest_path_best=os.path.join(dest_path_csv, "ensemble_"+formula+"_best.csv")

            my_df.to_csv(dest_path_all,index=False)

            my_df=my_df.drop_duplicates(subset=["sym"])
            my_df.to_csv(dest_path_best,index=False)
        else:
            for model_name in Model_name_list:
                my_df=Tools.data_reduction(path,model_name)
                dest_path_csv=os.path.join(path,"Calc_report",model_name)
                dest_path_all=os.path.join(dest_path_csv, model_name+"_"+formula+"_all_reduced.csv")
                dest_path_best=os.path.join(dest_path_csv, model_name+"_"+formula+"_best_reduced.csv")

                my_df.to_csv(dest_path_all,index=False)

                my_df=my_df.drop_duplicates(subset=["sym"])
                my_df.to_csv(dest_path_best,index=False)

    if(CheckNew==True):
        DBsearch.find_unique_data(formula,path)
        if(top_n!=None):
            tag="all_OnlyNew.csv"
            DBsearch.copy_best_data(path,tag,top_n)
        elif(top_c!=None):
            tag="all.csv"
            DBsearch.copy_best_data(path,tag,top_c)
        else:
            pass
    else:
        if(top_n!=None):
            tag="all_OnlyNew.csv"
            DBsearch.copy_best_data(path,tag,top_n)
        elif(top_c!=None):
            tag="all.csv"
            DBsearch.copy_best_data(path,tag,top_c)
        else:
            pass

    if(StoreData==True):
        from PNcsp_Plus.db import DBconnector
        if calculator==None:
            print("Error: Calculator is missing! Enter the calculator (MACE, M3GNet, ALIGNN).")
            return(None)
        if dim==None:
            print("Error: Dimension is missing! Enter the dimension.")
            return(None)
        # if project==None:
        #     print("Error: Project is missing! Enter the project.")
        #     return(None)
        DBconnector.upload_data(data_path, formula, calculator, dim, project)

def main(args_list=None):

    def parse_None_float(value):
        if value==None:
            return None
        if type(value) == str:
            if value.lower() == 'none':
                return None
            try:
                # Convert valid numeric strings to float
                return float(value)
            except ValueError:
                raise argparse.ArgumentTypeError(f"Invalid filter value: '{value}'. Expected a number or 'none'.")
    
    def parse_None_all_int(value):
        if value==None:
            return None
        if type(value) == str:
            if value.lower() == 'none':
                return None
            if value.lower() == 'all':
                return "all"
            try:
                # Convert valid numeric strings to float
                return int(value)
            except ValueError:
                raise argparse.ArgumentTypeError(f"Invalid filter value: '{value}'. Expected an int, 'all' or 'none'.")

    def parse_None_string(value):
        if value==None:
            return None
        if type(value) == str:
            if value.lower() == 'none':
                return None
            else:
                return value
        
    parser = argparse.ArgumentParser(prog="PNcsp",description= "PNcsp: A PN similarty based initial structure generator.")
    parser.add_argument('formula')
    parser.add_argument('-n','--neighbor',default=1,help="Order of neighbors to be considered in the similarity search. (default: 1)")
    parser.add_argument('-f','--filter',default=0,help="Selected neighbors are limited to those below the energy filter value. (default: 0) unit: [eV/atom]. Use \"None\" for no filter.") 
    parser.add_argument('-ts','--time_sleep',default=None,help="Set sleep time between queries. Excessive number of queries may cause the server to halt.(default: \"None\")")
    parser.add_argument('-on','--online',default="False",help="Sets online (True) or offline (False) search in OQMD (default: False). For offline seach, you should download and set up offline OQMD database. See https://oqmd.org/download/.")
    parser.add_argument('-calc','--calculator',default=None,help="Sets calculator [M3GNet, ALIGNN, MACE, ensemble]. Calculators are applied to all available neighbors. (default: None).")
    parser.add_argument('-out','--output_dir',default=".",help="Sets output directory. Enter full path. (default: current directory).")
    parser.add_argument('--BlockSearch',help="Blocks search. In case you want to use only calculator but not search feature, use this flag.",action='store_true')
    parser.add_argument('--Relax',help="Sets Structure relaxation before ML evaluation.",action='store_true')
    parser.add_argument('-db','--database',default="OQMD",help="Sets data source [OQMD, MP, MPDS]. (default: OQMD).")
    parser.add_argument('--CheckNew',help="Check if found structures have been already reported in OQMD and MP.",action='store_true')
    parser.add_argument('-top_n','--top_n_new',default=None,help="Copies the top-n evaluated new structures, ranked by the GNN evaluation, to the Best_Structures folder if available [int, \"all\"]. (default: None). Use this option together with --CheckNew, or in a subsequent run after a run performed with --CheckNew. ")
    parser.add_argument('-top_c','--top_n_calc',default=None,help="Copies the top-n evaluated structures, ranked by the GNN evaluation, to the Best_Structures folder if available [int, \"all\"]. (default: None). Use this option together with --CheckNew, or in a subsequent run after a run performed with --CheckNew. ")
    parser.add_argument('-project','--project',default=None,help="Assing project name. This option is only used together with --StoreData.")
    parser.add_argument('-dim','--dimension',default=None,help="Assing dimension of the crystalline system. This option is only used together with --StoreData. (default: None).")
    parser.add_argument('--ReduceData',help="Removes duplicates.",action='store_true')
    parser.add_argument('--Reverse',help="Searces neigborhood of given formula",action='store_true')
    parser.add_argument('--StoreData',help="Stores generated data in local DB. For this option, you should have a local DB. ",action='store_true')

    args = parser.parse_args(args_list)

    formula=args.formula
    N_neig=int(args.neighbor)
    online = False if args.online == "False" else True
    calculator=parse_None_string(args.calculator)
    database=args.database
    BlockSearch=args.BlockSearch
    ReduceData=args.ReduceData
    CheckNew=args.CheckNew
    Reverse=args.Reverse
    Relax=args.Relax
    StoreData=args.StoreData
    project=args.project
    dim=args.dimension

    if dim:
        dim=int(dim)
    
    time_sleep=parse_None_float(args.time_sleep)
    E_filter=parse_None_float(args.filter)

    if(args.output_dir[-1]=="/"):
        data_path=args.output_dir[:-1]
    else:
        data_path=args.output_dir

    top_n=parse_None_all_int(args.top_n_new)
    top_c=parse_None_all_int(args.top_n_calc)

    run(formula,N_neig=N_neig,E_filter=E_filter,time_sleep=time_sleep,online=online,calculator=calculator,data_path=data_path,
        BlockSearch=BlockSearch,Relax=Relax,database=database,CheckNew=CheckNew,top_n=top_n,
        top_c=top_c,project=project,dim=dim,ReduceData=ReduceData,Reverse=Reverse,StoreData=StoreData)

if __name__=='__main__':
    main()
