def offline_order(Comp_list):
    import re
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

def get_data_OQMD(Comp_list,neigh_list,Energy_filter,timer):
    import qmpy

    All_list=[]
    Comp_list_ordered=offline_order(Comp_list)
    print(Comp_list)
    a=[]
    for i in range(len(Comp_list_ordered)):
        raw_data = qmpy.Entry.objects.filter(composition_id=Comp_list_ordered[i])
        if(len(raw_data)==0):
            print(Comp_list_ordered[i],"--> no structure")
            continue

        list_of_data={"data":[]}
        for k in range(len(raw_data)):
            comp=raw_data[k]
            delta_e=comp.energy
            if delta_e<=Energy_filter:
                list_of_data["data"].append({"name":comp.name, "spacegroup":comp.spacegroup, "cell":comp.structure.cell, "sites":comp.structure.atoms, "Neigh":neigh_list[i],"Original":Comp_list[i]})
        if list_of_data["data"]!=[]:
            All_list.append(list_of_data)
            print(Comp_list_ordered[i])
        else:
            print(Comp_list_ordered[i],"--> no structure")
    if(All_list==[]):
        print("Warning: No candidates were found! TERMINATED!")
        exit(0)
    return All_list

def create_prototype_OQMD(All_list,exchange_dict,formula):
    import os
    from ase.spacegroup import crystal
    from ase.io import read, write

    path="./output_"+formula+"/"
    for num1 in range(len(All_list)):

        neigh=All_list[num1]['data'][0]['Neigh']
        dest_path=path+str(neigh)+"_Neigh/"

        compound=All_list[num1]
        for num2 in range(len(All_list[num1]['data'])):
            name=compound['data'][num2]['name']
            spacegroup=compound['data'][num2]['spacegroup'].symbol
            cell=compound['data'][num2]['cell']
            sites=compound['data'][num2]['sites']

            site_list=[]
            elem_list=[]

            for i in range(len(sites)):
                site=tuple(str(sites[i]).split(' @ ')[1].split(' '))
                site_list.append(site)
                elem=str(sites[i]).split(' @ ')[0]
                elem_list.append(elem)

            for i in range(len(elem_list)):
                elem_list[i]=exchange_dict[elem_list[i]]

            struct = crystal(elem_list, site_list, cell=cell,size=(1,1,1))

            if not os.path.exists(dest_path):
                os.makedirs(dest_path)
            
            write(dest_path+formula+"_"+name+"_"+str(spacegroup).replace("/","")+"_"+str(num2)+'.cif',struct)
