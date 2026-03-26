def get_data_MP(comp_list,Energy_filter):
    from mp_api.client import MPRester

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
        All_list=[]
        for comp in comp_list:
            data=mpr.summary.search(formula=comp,fields=fields)

            data_neg = [dat for dat in data if dat.formation_energy_per_atom <= Energy_filter]
            if data_neg==[]:
                continue
            All_list.append(data_neg)
        return All_list
    
def create_prototype_MP(All_list, exchange_dict, formula, neigh,data_path):
    import os
    from itertools import product
    from ase.spacegroup import crystal
    import re
    from math import gcd
    from functools import reduce
    import shutil
    from ase.io import write
    from pymatgen.core import Structure

    def write_CIF(elem_list_replaced, dest_path, formula, name, spacegroup, num2, ind_ext):
        struct = crystal(elem_list_replaced, site_list, cell=cell, size=(1, 1, 1))
        filename = f"{dest_path}{formula}_{name}_{str(spacegroup).replace('/','')}_{num2}_{ind_ext}.cif"
        write(filename, struct)
    
    def write_CIF2(struct, dest_path, formula, name, spacegroup, num2, ind_ext):

        filename = f"{dest_path}{formula}_{name}_{str(spacegroup).replace('/','')}_{num2}_{ind_ext}.cif"
        # write(filename, struct)
        struct.to(fmt="cif", filename=filename)

    def reduce_list(nums):
        g = reduce(gcd, nums)
        return [n // g for n in nums]

    
    def parse_formula_counts(formula):
        matches = re.findall(r'[A-Z][a-z]*|\d+', formula)
        return {matches[i]: int(matches[i+1]) for i in range(0, len(matches), 2)}

    path = data_path+"/output_"+formula+"/"
    expected_counts = parse_formula_counts(formula)
    
    dest_path = f"{path}{neigh}_Neigh/"
    
    if not os.path.exists(dest_path):
        os.makedirs(dest_path)
    else:
        shutil.rmtree(dest_path)
        os.makedirs(dest_path)

    for num1, compound in enumerate(All_list):
        for num2, entry in enumerate(compound):
            name = entry.formula_pretty
            structure = entry.structure
            spacegroup = entry.symmetry.symbol
            site_list=[site.frac_coords for site in entry.structure]
            cell = entry.structure.lattice.parameters

            # Collect element list
            elem_list=[site.specie.symbol for site in structure]
            elem_list_org=[site.specie.symbol for site in structure]

            # Identify ambiguous elements (those with >1 possible substitution)
            ambiguous_elems = list({e for e in elem_list if len(exchange_dict[e]) > 1})

            ind_ext = 0
            # # print(exchange_dict)
            if not ambiguous_elems:
                # No ambiguity → just replace with the unique option
                elem_list_new = [exchange_dict[e][0] for e in elem_list]
                for e in elem_list:
                    structure.replace_species({e:exchange_dict[e][0]})   
                # print(name, ambiguous_elems, "\n", elem_list_new, "\n", elem_list)
                write_CIF2(structure, dest_path, formula, name, spacegroup, num2, ind_ext)
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
            #         # print(name, ambiguous_elems, "\n", elem_list_new, "\n", elem_list)
                    continue

                
                unique_elems = list(set(elem_list_new))
                counts = reduce_list([elem_list_new.count(e) for e in unique_elems])

                match = True
                for e, cnt in zip(unique_elems, counts):
                    if expected_counts.get(e, -1) != cnt:
                        match = False
                        break

                if match:
                    write_CIF(elem_list_new, dest_path, formula, name, spacegroup, num2, ind_ext)

                    ind_ext += 1