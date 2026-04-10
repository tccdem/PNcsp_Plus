#!/bin/bash

# Path to your PNcsp Python script or executable
OUT="test_data"

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
    # echo "***** $compound neigh=3 ******";
    # "${PNCSP_CMD[@]}" -n 3;
    # echo "***** $compound neigh=4 ******";
    # "${PNCSP_CMD[@]}" -n 4;
    echo "***********************";
    echo "----------------------------";
done
echo "+++ date 2 +++";
date

for compound in "${COMPOUNDS[@]}"; do
    PNCSP_CMD=(python PNcsp.py "$compound" -out "$OUT")
    echo "***** $compound ******";
    "${PNCSP_CMD[@]}" -calc MACE --BlockSearch;
    echo "***********************";
    echo "----------------------------";
done

echo "+++ date 3 +++";
date

for compound in "${COMPOUNDS[@]}"; do
    PNCSP_CMD=(python PNcsp.py "$compound" -out "$OUT")
    echo "***** $compound ******";
    "${PNCSP_CMD[@]}" --CheckNew;
    echo "***********************";
    echo "----------------------------";
done
