#!/bin/bash

# Path to your PNcsp Python script or executable
OUT="dev/test_data"

# List of compounds
COMPOUNDS=(
Ho1Pb1
)

echo "+++ date 1 +++";
date 
# Neighbor search + Generation of prototypes
echo "***** Neighbor Search ******";
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
    echo "----------------------------";
done
echo "+++ date 2+++";
date 

# GNN evaluation + Data reduction (internal) (generate csv files under Calc_report)
echo "***** GNN evaluation + Data reduction ******";
for compound in "${COMPOUNDS[@]}"; do
    PNCSP_CMD=(python PNcsp.py "$compound" -out "$OUT")
    echo "***** $compound ******";
    "${PNCSP_CMD[@]}" --BlockSearch -calc MACE --ReduceData;
    echo "----------------------------";
done

echo "+++ date 3+++";
date 

# Detection of new structures (external) + Copy Data 
echo "***** Detection of new structures + Copy Data ******";
for compound in "${COMPOUNDS[@]}"; do
    PNCSP_CMD=(python PNcsp.py "$compound" -out "$OUT")
    echo "***** $compound ******";
    # "${PNCSP_CMD[@]}" --BlockSearch -top_c all;
    "${PNCSP_CMD[@]}" --BlockSearch --CheckNew -top_c all;
    echo "***********************";
    echo "----------------------------";
done

# # MACE relaxation to copied data
# echo "***** MACE Relaxation ******";
# for compound in "${COMPOUNDS[@]}"; do
#     RELAXER_CMD=(python MACE_relaxer.py "$compound" -out "$OUT")
#     "${RELAXER_CMD[@]}";
#     echo "***********************";
# done

# # Generation of CASTEP inputs
# echo "***** Generation of CASTEP inputs ******";
# for compound in "${COMPOUNDS[@]}"; do
#     echo $compound;
#     CIF_DIR="$OUT/output_$compound/Calc_report/MACE/CheckNew_report/Best_Structures/opt"

#     ENTERPRISE_CMD=(python dev/enterprise/castep-automation-wconfig.py --cif_dir "$CIF_DIR")
#     "${ENTERPRISE_CMD[@]}";

#     # ENTERPRISE_CMD=(./autocasp --cif_dir "$CIF_DIR")
#     # "${ENTERPRISE_CMD[@]}";
#     echo "***********************";
#     echo "----------------------------";
# done
