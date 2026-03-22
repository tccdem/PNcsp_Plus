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

def get_data_OQMD(Comp_list,neigh_list,Energy_filter):
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

def create_prototype_OQMD(All_list, exchange_dict, formula,data_path):
    import os
    from ase.spacegroup import crystal
    from ase.io import write
    from itertools import product
    import re
    from math import gcd
    from functools import reduce
    import shutil

    def write_CIF(elem_list_replaced, dest_path, formula, name, spacegroup, num2, ind_ext):
        struct = crystal(elem_list_replaced, site_list, cell=cell, size=(1, 1, 1))
        
        filename = f"{dest_path}{formula}_{name}_{str(spacegroup).replace('/','')}_{num2}_{ind_ext}.cif"
        write(filename, struct)
    
    def reduce_list(nums):
        g = reduce(gcd, nums)
        return [n // g for n in nums]

    
    def parse_formula_counts(formula):
        matches = re.findall(r'[A-Z][a-z]*|\d+', formula)
        return {matches[i]: int(matches[i+1]) for i in range(0, len(matches), 2)}

    path = data_path+"/output_"+formula+"/"
    expected_counts = parse_formula_counts(formula)

    neigh = All_list[0]['data'][0]['Neigh']
    dest_path = f"{path}{neigh}_Neigh/"

    if not os.path.exists(dest_path):
            os.makedirs(dest_path)
    else:
        shutil.rmtree(dest_path)
        os.makedirs(dest_path)

    for num1, compound in enumerate(All_list):
        for num2, entry in enumerate(compound['data']):
            name = entry['name']
            spacegroup = entry['spacegroup'].symbol
            cell = entry['cell']
            sites = entry['sites']

            # Collect element list and site list
            site_list, elem_list, elem_list_org = [], [], []
            for s in sites:
                elem, coords = str(s).split(' @ ')
                site_list.append(tuple(coords.split()))
                elem_list.append(elem)
                elem_list_org.append(elem)

            # Identify ambiguous elements (those with >1 possible substitution)
            ambiguous_elems = list({e for e in elem_list if len(exchange_dict[e]) > 1})

            ind_ext = 0
            # print(exchange_dict)
            if not ambiguous_elems:
                # No ambiguity → just replace with the unique option
                elem_list_new = [exchange_dict[e][0] for e in elem_list]
                # print(name, ambiguous_elems, "\n", elem_list_new, "\n", elem_list)
                write_CIF(elem_list_new, dest_path, formula, name, spacegroup, num2, ind_ext)
                ind_ext += 1
                continue

            # Build all possible substitution combinations for the ambiguous_elems elements
            replacement_options = [exchange_dict[e] for e in ambiguous_elems]

            for combo in product(*replacement_options):
                # Map chosen replacements to the ambiguous elements
                replacement_map = dict(zip(ambiguous_elems, combo))
                
                elem_list_new = [replacement_map.get(e, e) for e in elem_list]
                # Force unique mappings (non-ambiguous elements get fixed substitution)
                for i, e in enumerate(elem_list_new):
                    if len(exchange_dict[e]) == 1:
                        elem_list_new[i] = exchange_dict[e][0]

                # Skip if substitution reduces unique element count
                if len(set(elem_list_new)) < len(set(elem_list_org)):
                    # print(name, ambiguous_elems, "\n", elem_list_new, "\n", elem_list)
                    continue

                
                unique_elems = list(set(elem_list_new))
                counts = reduce_list([elem_list_new.count(e) for e in unique_elems])

                match = True
                for e, cnt in zip(unique_elems, counts):
                    if expected_counts.get(e, -1) != cnt:
                        match = False
                        break

                if match:
                    # print(name, ambiguous_elems, "\n", elem_list_new, "\n", elem_list)
                    write_CIF(elem_list_new, dest_path, formula, name, spacegroup, num2, ind_ext)
                    ind_ext += 1