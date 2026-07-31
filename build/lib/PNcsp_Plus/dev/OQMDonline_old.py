# OQMD
def get_data_OQMD(Comp_list,neigh_list,Energy_filter,timer):
    import qmpy_rester as qr
    import time

    All_list=[]
    for i in range(len(Comp_list)):
        with qr.QMPYRester() as q:
            if(Energy_filter=="none"):
                kwargs = {
                    'composition': {Comp_list[i]},
                    'format': 'json',
                    }
            else:
                kwargs = {
                    'composition': {Comp_list[i]},
                    'format': 'json',
                    'delta_e': "<"+str(Energy_filter),
                    }
            list_of_data = q.get_oqmd_phases(False,**kwargs, verify=False)

            if list_of_data is None:
                print("\nWARNING!")
                print("--------")
                print("Time exceed ! Wait for a while and use timer with higher value !!!")
                break
            
            if(list_of_data['data']==[]):
                print(Comp_list[i],"--> no structure")
                if(timer!="none"):
                    time.sleep(timer)
                continue
            
            for ind in range(len(list_of_data['data'])):
                list_of_data['data'][ind]['Neigh']=neigh_list[i]
                # list_of_data['data'][ind]['Original']=''.join([i for i in Comp_list[i] if not i.isdigit()])
                # list_of_data['data'][ind]['Original']=Comp_list[i].replace("1","")
                list_of_data['data'][ind]['Original']=Comp_list[i]


            All_list.append(list_of_data)
            print(Comp_list[i])
        if(timer!="none"):
            time.sleep(timer)
    if(All_list==[]):
        print("Warning: No candidates were found! TERMINATED!")
        exit(0)
    return All_list

def create_prototype_OQMD(All_list,exchange_dict,formula,data_path):
    import os
    from ase.spacegroup import crystal
    from ase.io import read, write

    path=data_path+"/output_"+formula+"/"
    for num1 in range(len(All_list)):

        neigh=All_list[num1]['data'][0]['Neigh']
        dest_path=path+str(neigh)+"_Neigh/"

        compound=All_list[num1]
        for num2 in range(len(All_list[num1]['data'])):
            name=compound['data'][num2]['name']
            # elems=re.findall('[A-Z][^A-Z]*', name)
            spacegroup=compound['data'][num2]['spacegroup']
            unit_cell=compound['data'][num2]['unit_cell']

            sites=compound['data'][num2]['sites']
            site_list=[]
            elem_list=[]
            for i in range(len(sites)):
                site=tuple(sites[i].split(' @ ')[1].split(' '))
                site_list.append(site)
                elem=sites[i].split(' @ ')[0]
                elem_list.append(elem)

            for i in range(len(elem_list)):
                elem_list[i]=exchange_dict[elem_list[i]]

            # # IF exchange_dict[elem_list[i]] values are the same, detect it and do related operations.
            # list_dup=[]
            # for i in range(len(elem_list)-1):
            #     dup=[]
            #     for j in range(i+1,len(elem_list)):
            #         if(elem_list[i]==elem_list[j]):
            #             dup.append(j)
            #     if(dup!=[]):
            #         dup.append(i)
            #         list_dup.append(dup)

            # not_included = []
            # for element in list(exchange_dict.values()):
            #     if element not in elem_list:
            #         not_included.append(element)
            # print("not_included:",not_included)
            # print("elem_list:",elem_list)

            struct = crystal(elem_list, site_list, cell=unit_cell,size=(1,1,1))
            if not os.path.exists(dest_path):
                os.makedirs(dest_path)
            
            write(dest_path+formula+"_"+name+"_"+spacegroup.replace("/","")+"_"+str(num2)+'.cif',struct)
