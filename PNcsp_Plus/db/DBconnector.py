from pymatgen.core import Structure
from tccdem_api import MaterialDatabaseAPI
import pandas as pd
import os


def get_path(dest_path, formula,program):
    struc_path= os.path.join(dest_path,"output_"+formula+"/Calc_report/"+program+"/CheckNew_report/Best_Structures/")
    info_path= os.path.join(dest_path,"output_"+formula+"/Calc_report/"+program+"/CheckNew_report/"+program+"_"+formula+"_all.csv")
    
    if not ((os.path.exists(struc_path)) or (os.path.exists(info_path))):
        print("Error: Calculations are missing! Complete to store the data!!")
        # exit(0)
        return None, None
    """
    struc_path= "/home/cem/PNcsp_Plus/PNcsp_Plus/dev/test_data/output_Ho1Pb1/Calc_report/MACE/CheckNew_report/Best_Structures/"
    info_table= "/home/cem/PNcsp_Plus/PNcsp_Plus/dev/test_data/output_Ho1Pb1/Calc_report/MACE/CheckNew_report/MACE_Ho1Pb1_all.csv"
    """
    return struc_path, info_path

def upload_data(dest_path,fomula,program,dim,project):
    # Initialize your DB API client
    db = MaterialDatabaseAPI()

    struc_path,info_path=get_path(dest_path,fomula,program)

    if not struc_path:
        return None

    files=os.listdir(struc_path+"/")
    for file in files:
        structure = Structure.from_file(struc_path+file)

        df = pd.read_csv(info_path)
        file_name = file.replace(".cif","")[:-4]
        print(file_name)
        df = df[df['CIF_Name'].str.contains(file_name)]
        delta_e = df.Energy.item()
        is_struct_new = df.is_struc_new.item()
        is_sg_new = df.is_sym_new.item()
        program = info_path.split("/")[-3]
        properties_payload = {
            "delta_e": {
                "value": delta_e,
                "unit": "eV/atom",
                "program": program
            },
        }

        # # 3. Fire the pipeline directly into the database
        db.upload_data(
            project=project,
            structure_obj=structure,
            properties_dict=properties_payload,
            dimension=dim,
            prototype="not_known",
            is_struct_new=is_struct_new,        # Explicitly marked as new prediction
            is_sg_new=is_sg_new,
            relaxation=None
        )
    print("Transfer to the TCCDEM_DB is succesfully done!!")