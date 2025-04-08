import os
import glob
import csv
from qpath_equilibrator_pathway_function import *

folder_path = "./result(20)/output(ACCOA)1"
csv_output = './result(20)//output(ACCOA)1/Pathway_thermodynamics_summary_results.csv' 

txt_files = glob.glob(os.path.join(folder_path, "*.txt"))

results = [["File", "deltaG", "MDF"]]  

for txt_file in txt_files:
    tsv_file = txt_file.replace(".txt", ".tsv")

    with open(txt_file, 'r', encoding='utf-8') as f_in:
        lines = [line.strip().split('\t')[:3] for line in f_in]  

    header = ["Reaction ID", "Fluxes", "Equation(ID)"]
    lines = [header] + lines

    with open(tsv_file, 'w', encoding='utf-8', newline='') as f_out:
        for line in lines:
            f_out.write('\t'.join(line) + '\n')

    substrate = 'CARBON__45__DIOXIDE_c'
    delmetcompartment = True
    outfile = tsv_file.replace(".tsv", "-sbtab.tsv") 
    pathway_transfer_sbtab(substrate, tsv_file, outfile, delmetcompartment)

    alldgfile = tsv_file.replace(".tsv", "-alldg.tsv") 
    mdf_figurefile = tsv_file.replace(".tsv", "-mdf.png")
    sta_figurefile = tsv_file.replace(".tsv", "-summdg.png")

    deltaG, MDF = cal_pathway_dg(outfile, sta_figurefile, mdf_figurefile, alldgfile)

    results.append([os.path.basename(txt_file), deltaG, MDF])
    results[1:] = sorted(results[1:], key=lambda x: x[0])  

with open(csv_output, 'w', encoding='utf-8', newline='') as f_out:
    writer = csv.writer(f_out)
    writer.writerows(results)