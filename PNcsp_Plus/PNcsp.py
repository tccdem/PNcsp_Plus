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

from calculator import ML


def blockPrint():
    sys.stdout = open(os.devnull, 'w')
def enablePrint():
    sys.stdout = sys.__stdout__

def get_Symbol(PN):
    with open("./db/data/Z_PN_Elem_extended_MD.csv",'r') as file:
        Comp_list={}
        for line in file:
            if 'PN' in line:
                continue
            info=line.strip().replace('\ufeff','').split(',')
            Comp_list[info[0]]=info[1]
    # print(Comp_list)
    return(Comp_list[str(PN)])

def get_PN(Symb):
    with open("./db/data/Z_PN_Elem_extended_MD.csv",'r') as file:
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
                # print("\nWARNING: PN distance between some of constituent elemens are lower than N_neig. This may cause peculiar or incomplete outputs.\n")
                print("\nWARNING: The neighbors of the ",elems[i],"and",elems[j] ,"overlap. This is not a problem, but it is advisable to exercise caution when examining the output CIFs.")
    print("\n")
    
    # # Exchange Look-up Table
    # exchange_dict={}
    # for i in range(len(PN_new_all)):
    #     for j in range(len(PN_new_all[i])):
    #         key=get_Symbol(PN_new_all[i][j])
    #         # print(key,":",elems[i])
    #         exchange_dict[key]=elems[i]
    # print("\nexchange_dict: ",exchange_dict)

    # # Exchange Look-up Table for PN
    # exchange_dict_PN={}
    # for i in range(len(PN_new_all)):
    #     for j in range(len(PN_new_all[i])):
    #         key=PN_new_all[i][j]
    #         # print(key,":",elems[i])
    #         exchange_dict_PN[key]=PNs[i]
    # print("\nexchange_dict_PN: ",exchange_dict_PN)

    # New Exchange Look-up Table
    exchange_dict={}
    for i in range(len(PN_new_all)):
        for j in range(len(PN_new_all[i])):
            key=get_Symbol(PN_new_all[i][j])
            # print(key,":",elems[i])
            if exchange_dict.get(key) is not None:
                exchange_dict[key].append(elems[i])
            else:
                exchange_dict[key]=[elems[i]]
    print("\nexchange_dict: ",exchange_dict)

    # Exchange Look-up Table for PN
    exchange_dict_PN={}
    for i in range(len(PN_new_all)):
        for j in range(len(PN_new_all[i])):
            key=PN_new_all[i][j]
            # print(key,":",elems[i])
            if exchange_dict_PN.get(key) is not None:
                exchange_dict_PN[key].append(PNs[i])
            else:
                exchange_dict_PN[key]=[PNs[i]]
    # print("\nexchange_dict_PN: ",exchange_dict_PN)


    res=np.array(list(itertools.product(*PN_new_all)))
    
    # # Drop original formula
    # for i in range(len(res)):
    #     if(np.array_equal(res[i], PNs)):
    #         res=np.concatenate((res[:i],res[i+1:]), axis=0)
    #         break

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

    # print("\nOrder_list:", dist_list_ordered)
    # print("\nPN  combinations:\n",res_ordered)

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

    # print("\nOrder_list:", dist_list_ordered)
    # print("\nPN  combinations and Symbol list")
    # for i in range(len(res_ordered)):
    #     print(res_ordered[i],Symbol_list[i])
    return Symbol_list,dist_list_ordered,exchange_dict

def categorize(N_neig,formula,data_path):
    import shutil  
    import os
    
    path=data_path+"/output_"+formula+"/"+str(N_neig)+"_Neigh/"
    cifs=[f for f in os.listdir(path) if ".cif" in f]
    print("Total Number of CIF files:",len(cifs))
    for cif in cifs:
        source_path=path+cif

        sym=cif.split("_")[2]
        dest_path=path+sym+"/"

        if not os.path.exists(dest_path):
            os.makedirs(dest_path)

        dest = shutil.move(source_path, dest_path)  

def show_config(formula,N_neig,E_filter,timer,online,calculator,database,BlockSearch,Relaxer,data_path):
    print("\nProgram Configuration")
    print("---------------------")
    print("Query formula:\t",formula,"\nNeighbor order:\t",N_neig,
          "\nEnergy filter:\t",E_filter,"\nSleep timer:\t",timer,
          "\nOnline: \t",online,"\nCalculator: \t",calculator,"\Data source: \t",database,
          "\nBlockSearch: \t",BlockSearch,"\nRelaxer: \t",Relaxer,
          "\nOutputDir: \t",data_path)
    print("---------------------\n")

def main():
    parser = argparse.ArgumentParser(prog="PNcsp",description= "PNcsp: A PN similarty based initial structure generator.")
    parser.add_argument('formula')
    parser.add_argument('-n','--neighbor',default=1,help="Order of neighbors to be considered in the similarity search. (default: 1)")
    parser.add_argument('-f','--filter',default=0,help="Selected neighbors are limited to those below the energy filter value. (default: 0) unit: [eV/atom]. Use \"none\" for no filter.") 
    parser.add_argument('-t','--time_sleep',default="none",help="Set sleep time between queries. Excessive number of queries may cause the server to halt.(default: \"none\")")
    parser.add_argument('-o','--online',default="False",help="Sets online (True) or offline (False) search in OQMD (default: False). For offline seach, you should download and set up offline OQMD database. See https://oqmd.org/download/.")
    parser.add_argument('-calc','--calculator',default="None",help="Sets calculator [M3GNet, ALIGNN, MACE, ensemble]. Calculators are applied to all available neighbors. (default: None).")
    parser.add_argument('-out','--output_dir',default=".",help="Sets output directory. Enter full path. (default: current directory).")
    parser.add_argument('--BlockSearch',help="Blocks search. In case you want to use only calculator but not search feature, use this flag.",action='store_true')
    parser.add_argument('--Relax',help="Sets Structure relaxation before ML evaluation.",action='store_true')
    parser.add_argument('-db','--database',default="OQMD",help="Sets data source [OQMD, MP, MPDS]. (default: OQMD).")

    args = parser.parse_args()

    formula=args.formula

    N_neig=int(args.neighbor)
    online = False if args.online == "False" else True
    calculator=args.calculator
    database=args.database
    BlockSearch=args.BlockSearch
    Relax=args.Relax
    
    if(args.time_sleep =="none"):
        time_sleep=args.time_sleep
    else:
        time_sleep=int(args.time_sleep)

    if(args.filter =="none"):
        E_filter=args.filter
    else:
        E_filter=int(args.filter)

    if(args.output_dir[-1]=="/"):
        data_path=args.output_dir[:-1]
    else:
        data_path=args.output_dir

    show_config(formula=formula,N_neig=N_neig,E_filter=E_filter,timer=time_sleep,online=online,calculator=calculator,database=database,BlockSearch=BlockSearch,Relaxer=Relax,data_path=data_path)
    
    if(BlockSearch!=True):
        res,neigh_list,exchange_dict=get_Neig(formula=formula,N_neig=N_neig)

        if(database=="OQMD"):
            if(online==True):
                from db import OQMDonline
                All_list=OQMDonline.get_data_OQMD(res,neigh_list,Energy_filter=E_filter,timer=time_sleep)
                OQMDonline.create_prototype_OQMD(All_list,exchange_dict,formula=formula,data_path=data_path)
            else:
                from db import OQMDoffline
                All_list=OQMDoffline.get_data_OQMD(res,neigh_list,Energy_filter=E_filter)
                OQMDoffline.create_prototype_OQMD(All_list,exchange_dict,formula=formula,data_path=data_path)
        elif(database=="MP"):
            from db import MPonline
            All_list=MPonline.get_data_MP(res,Energy_filter=E_filter)
            MPonline.create_prototype_MP(All_list, exchange_dict, formula=formula, neigh=N_neig,data_path=data_path)
        elif(database=="MPDS"):
            print("WARNING: MPDS implementation is under construction. Choose another data source and run again.")
            exit(0)
            
        print("TERMINATED SUCCESFULLY!\n")
        categorize(N_neig=N_neig,formula=formula,data_path=data_path)

    path=data_path+"/output_"+formula+"/"
    if(calculator=="M3GNet"):
        # ML.M3GNet_calc("./","./output_"+formula+"/Calc_report/")
        path_results=data_path+"/output_"+formula+"/Calc_report/M3GNet/"
        if(Relax==True):
            ML.M3GNet_calc(formula,path,path_results,relax=True)
        else:
            ML.M3GNet_calc(formula,path,path_results,relax=False)
    elif(calculator=="ALIGNN"):
        path_results=data_path+"/output_"+formula+"/Calc_report/ALIGNN/"
        if(Relax==True):
            ML.ALIGNN_calc(formula,path,path_results,relax=True)
        else:
            ML.ALIGNN_calc(formula,path,path_results,relax=False)
    elif(calculator=="MACE"):
        path_results=data_path+"/output_"+formula+"/Calc_report/MACE/"
        if(Relax==True):
            ML.MACE_calc(formula,path,path_results,relax=True)
        else:
            ML.MACE_calc(formula,path,path_results,relax=False)
    elif(calculator=="ensemble"):
        # path_results="./DATA/DATA_180_4NN/output_"+formula+"/Calc_report/ensemble/"
        path_results=data_path+"/output_"+formula+"/Calc_report/ensemble/"
        if(Relax==True):
            print("\nERROR: You can not run Majority Vote with relaxation.\n")
        else:
            # path_alignn="./DATA/DATA_180_4NN/output_"+formula+"/Calc_report/ALIGNN/"
            path_alignn=data_path+"/output_"+formula+"/Calc_report/ALIGNN/"
            ML.ALIGNN_calc(formula,path,path_alignn,relax=False)

            # path_mace="./DATA/DATA_180_4NN/output_"+formula+"/Calc_report/MACE/"
            path_mace=data_path+"/output_"+formula+"/Calc_report/MACE/"
            ML.MACE_calc(formula,path,path_mace,relax=False)

            # path_m3gnet="./DATA/DATA_180_4NN/output_"+formula+"/Calc_report/M3GNet/"
            path_m3gnet=data_path+"/output_"+formula+"/Calc_report/M3GNet/"
            ML.M3GNet_calc(formula,path,path_m3gnet,relax=False)
            
            N_model=2 # 2 Model: m3gnet, mace 3 Model: alignn, m3gnet, mace
            ML.ensemble_vote(formula,path_alignn,path_mace,path_m3gnet,path_results,N_model)
    else:
        pass

if __name__=='__main__':
    main()
