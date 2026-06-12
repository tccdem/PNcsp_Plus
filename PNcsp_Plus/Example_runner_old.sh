#!/bin/bash

# Path to your PNcsp Python script or executable
OUT="dev/test_data"

# List of compounds
COMPOUNDS=(
Al1Cu1
Mg1H1
)


echo "+++ date 1 +++";
date 
# Loop through compounds and run PNcsp
for compound in "${COMPOUNDS[@]}"; do
    PNCSP_CMD=(python PNcsp.py "$compound" -out "$OUT")

    echo "***** $compound neigh=1 ******";
    "${PNCSP_CMD[@]}" -n 1;
    echo "***** $compound neigh=2 ******";
    "${PNCSP_CMD[@]}" -n 2;
    echo "***** $compound neigh=3 ******";
    "${PNCSP_CMD[@]}" -n 3;
    echo "***** $compound neigh=4 ******";
    "${PNCSP_CMD[@]}" -n 4;
    echo "***********************";
    echo "----------------------------";
done
echo "+++ date 2 +++";
date

# for compound in "${COMPOUNDS[@]}"; do
#     PNCSP_CMD=(python PNcsp.py "$compound" -out "$OUT")
#     echo "***** $compound ******";
#     # "${PNCSP_CMD[@]}" -calc MACE --BlockSearch --CheckNew -top_n all;
#     "${PNCSP_CMD[@]}" --BlockSearch --CheckNew -top_n all;
#     echo "***********************";
#     echo "----------------------------";
# done

# echo "+++ date 3 +++";
# date

# echo "+++ date 3 +++";
# for compound in "${COMPOUNDS[@]}"; do
#     RELAXER_CMD=(python MACE_relaxer.py "$compound" -out "$OUT")
#     "${RELAXER_CMD[@]}";
#     echo "***********************";
# done
# echo "+++ date 4 +++";
# date

# for compound in "${COMPOUNDS[@]}"; do
#     echo $compound;
#     CIF_DIR="dev/test_data/output_Mg1H1/Calc_report/MACE/CheckNew_report/Best_Structures"

#     ENTERPRISE_CMD=(python dev/enterprise/castep-automation-wconfig.py --cif_dir "$CIF_DIR")
#     "${ENTERPRISE_CMD[@]}";

#     # ENTERPRISE_CMD=(./autocasp --cif_dir "$CIF_DIR")
#     # "${ENTERPRISE_CMD[@]}";
#     echo "***********************";
#     echo "----------------------------";
# done
