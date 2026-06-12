#!/bin/bash

# Path to your PNcsp Python script or executable
OUT="dev/CSPBench_data"

# List of compounds
COMPOUNDS=(
Na1Ga4
Sr1Ga4
Y1Mg1Cu1
La3Cu3Bi4
Cr1Fe1Co1Si1
Pu1Sn2
In1Fe1Co2
Sm1Fe5
U3Sb4Ir3
Cu2Si1Hg1S4
Ca1Cu3Pt4O12
In1Hg1
Yb3Co1
Ti1Ga3
Sm1H2
Fe2Cu6Sn1S8
La1Mn1Sb1O1
Ca1Ge2Pt1
Zr1Zn1Ni4
Pr1Cu3Ru4O12
Li1Mg2Ga1
Rb1Li3S2O9
Ce1Co1Si2
Ta2Fe1
Ce1Pb3
K1Br1F4
Th1B12
Co1Ni1Sn1
Zn1S1O4
Ca1Se1O3
Na2Cd1Pb1
Ba1Pr1Mn2O6
Ba1Nd2Co1O5
Ta1Ga1Co2
Ba1Cu1Te1F1
U2Mo1
K1Li6Ir1O6
Er1Ti2Ga4
Sn4Pd1
Ge3N4
Nd1Ga2Ni1
Al4Cu2O7
Sc1Ag1P2Se6
Ge4Rh1
Li1Mg1Sn1Pt1
Lu1In1Cu2
Yb1B12
Rb3Pr1Cl6
Ba2Yb1Nb1O6
K1Cu1Cl3
K2Na1Al1F6
Be1Pd2
Dy1Pb3
Th2Zn1
Ga2Cu1
K3Na1Fe1Cl6
K2Li1Cr1F6
Lu1Sn1Pd2
Ba4Na1B3N6
Yb1H3C1N3
Lu1Se1O3F1
Zr1Ga1
Sc1Cu1
Y1Ir1
Ba2Y1Ru1O6
Rb2Ti1O1F5
Tb1Cd2
Eu3Au2
Mn3Au1
Ga1Co1
Ho4Ga12Ni1
K1Er1S2
Mg3Au1
Th1Ga2
Ce1Al2B1Ru2
Ce1Nb1O4
Li1Mg1Sn1Au1
Ce1Cr2Si2C1
K1Rb2Sc1F6
Zr1Ta1N1O1
Sr1Ga1Cu2
Ba2U1Ni1O6
Dy1Cu1
Cu2Sn1Se3
Zn1Cd1Pt2
Hf1Co2Sn1
Ba2Tl2Cu1O6
Ho1Sn1Pt1
Tb5S7
Ce1Cu2Si2
La1F3
Sr1Ni1Sn3
Li1Ga1Si1
Co1Te1
Ba2Er1Sb1O6
Re1O3
Sm1Fe1As1O1
Ba2Pr1Sb1O6
La1Zn1Sb1O1
Cr3Ga1
K2Na1Ga1P2
Yb1Ga2
Be1Si1N2
Co4Ni1Sb12
Ga1P1
Y1Al3
Mg1In1Cu4
Ca1Ag1As1
Ba2Eu1Ta1O6
In1Pb2I5
La5Ag1Pb3
K2Li1In1F6
Sr2Mg1Ir1O6
Yb1Cu5
Ga2Os1
Dy1Pd1
Na1Pb2I1O6
Mg1Cu4Sn1
Lu1Ag1Pb1
Lu1Mn2Ge2
Er1Si2Au2
Ba2Bi1Sb1O6
Li2Ni1O2
Ho1Fe4Cu3O12
Cr1Te4Au1
Co3Sb4O6F6
Sr1Fe1As1F1
Li2Cu1Sn1
Zr5Al1Sb3
V1Cr1O4
Sr1Cu1S1F1
Ta4Ga1Te4Se4
Ba2Y1Ir1O6
Pr1Ni2B2C1
Cs1Cd1N3O6
Pr1C2
Sn1S1
Tm3Pt4
Ba2Sr1Te1O6
Dy1Cu2
Sm1Pd3S4
Pa1O1
Re1B4
Ba3Zn1N2O1
Hf2Ni1
Ce1Ga2
Er3Cu3Sb4
K2Na1In1F6
Sm1Ta1O4
Nb1P1Se1
Pr3Ir1
K2Na1Al1H6
La1B2Rh2C1
Ca1Cd2P2
Hf1Mn2
Mg2Be1B2Ir5
Y1Hg2
Nd1Ge2
Rb1P1H2O4
Mg1V4Sn1O12
Zr1Hg1
La1Cr4Cu3O12
Er1Co5
Y1Mn1O3
Ta1Be3
Ca3Sn1O1
Ti2Cd1
La1Cu1Te1O1
Nb3Si1
K3Mn1O4
Sr2Bi2Se3O2
Lu4Ga12Ni1
Mn1Bi1
Zn1C1O3
Ce1Mn4Cu3O12
K1As4I1O6
Ce2Sn1S5
Er1In1Cu2
Tb4Al1
Cs1U1F6
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
    "${PNCSP_CMD[@]}" --BlockSearch --CheckNew -top_c 5;
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

# Generation of CASTEP inputs
echo "***** Generation of CASTEP inputs ******";
for compound in "${COMPOUNDS[@]}"; do
    echo $compound;
    CIF_DIR="$OUT/output_$compound/Calc_report/MACE/CheckNew_report/Best_Structures/opt"

    ENTERPRISE_CMD=(python dev/enterprise/castep-automation-wconfig.py --cif_dir "$CIF_DIR")
    "${ENTERPRISE_CMD[@]}";

    # ENTERPRISE_CMD=(./autocasp --cif_dir "$CIF_DIR")
    # "${ENTERPRISE_CMD[@]}";
    echo "***********************";
    echo "----------------------------";
done
