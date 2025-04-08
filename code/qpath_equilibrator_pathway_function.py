from sbtab.SBtab import SBtabTable, SBtabDocument
import pandas as pd
import re
import os
from equilibrator_api import ComponentContribution, Q_
from equilibrator_pathway import ThermodynamicModel
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import pandas as pd
import json

def config():
    data = [
    ("algorithm", "MDF", "ECM, or MDF"),
    ("p_h", 7, ""),
    ("ionic_strength", "250 mM", ""),
    ("p_mg", 3, ""),
    ("dg_confidence", 0.95, ""),
    ]

    df = pd.DataFrame(data, columns=["!Option", "!Value", "!Comment"])

    config_table = SBtabTable.from_data_frame(
        df, table_id="Configuration", table_type="Config",table_name=None
    )
    return config_table

def replace_metabolites(reaction_formula,metacyc_df):
    # Divide the reaction equation into reactants and products according to "<=>"
    reactants, products = reaction_formula.split(" <=> ")

    # Separate the reactants and products into individual metabolites with a "+" sign
    reactants = reactants.split(" + ")
    products = products.split(" + ")

    # Traverse the reactants and replace the BioCyc number
    for i, reactant in enumerate(reactants):
        if reactant in list(metacyc_df.index):
            reactants[i] = metacyc_df.loc[reactant,'BioCyc']

    # Traverse the product and replace the BioCyc number
    for i, product in enumerate(products):
        if product in list(metacyc_df.index):
            products[i] = metacyc_df.loc[product,'BioCyc']

    # Reconstruct the replaced reaction equation
    new_reaction_formula = " + ".join(reactants) + " <=> " + " + ".join(products)

    return new_reaction_formula

def reaction(pathway_df,metacyc_df,delmetcompartment):
    # reaction
    reaction_data = []
    for ri in pathway_df.index:
        rxn_turple = ()
        equ = pathway_df.loc[ri,'Equation(ID)']

        # Convert the metabolite ID containing metacyc in the equation to biocyc ID
        equ = replace_metabolites(equ,metacyc_df)
        
        if delmetcompartment is True:
            equ = re.sub(r'(_[a-zA-Z]+)(?=\s|$)', '', equ) # Remove the compartments in the reaction _c,_p
        if delmetcompartment is False:
            equ = equ
            
        # Judging the reaction flux
        flux = pathway_df.loc[ri,'Fluxes']

        if flux > 0:
            rxn_turple = rxn_turple + (ri,equ)
            reaction_data.append(rxn_turple)
        else:
            equ = equ.split(' <=> ')[1] + ' <=> ' + equ.split(' <=> ')[0] # The flux is negative, and the reactant and product exchange positions in the reaction equation
            rxn_turple = rxn_turple + (ri,equ)
            reaction_data.append(rxn_turple)

    df = pd.DataFrame(reaction_data, columns=["!ID", "!ReactionFormula"])

    reaction_table = SBtabTable.from_data_frame(
        df, table_id="Reaction", table_type="Reaction"
    )
    return reaction_table

def compound(pathway_df,metacyc_df,kegg_df,delmetcompartment):
    # compound
    pathway_compound_all = []
    for ri in pathway_df.index:
        equ = pathway_df.loc[ri,'Equation(ID)']
        # Convert the metabolite ID containing metacyc in the equation to biocyc ID
        equ = replace_metabolites(equ,metacyc_df)
        
        if delmetcompartment is True:
            equ = re.sub(r'(_[a-zA-Z]+)(?=\s|$)', '',equ) # Remove reaction compartment _c,_p
        if delmetcompartment is False:
            equ = equ

        for mi in equ.split(' <=> ')[0].split(' + '):
            pathway_compound_all.append(mi.split(' ')[-1])
        for mi in equ.split(' <=> ')[1].split(' + '):
            pathway_compound_all.append(mi.split(' ')[-1])
        # rxn_met_list = re.findall(r'\w+__\w+|\w+', rxn)
    pathway_compound_all = list(set(pathway_compound_all))
    # print(len(pathway_compound_all))  
    # print(pathway_compound_all)

    compound_data = []
    for mi in pathway_compound_all:
        met_turple = ()

        if mi in list(metacyc_df['BioCyc']):
            met_turple = met_turple + (mi,'metacyc.compound:'+mi)
        elif mi in list(kegg_df.index):
            met_turple = met_turple + (mi,'kegg:'+mi)
        else:
            met_turple = met_turple + (mi,'bigg.metabolite:'+mi)
        compound_data.append(met_turple)
    # print(compound_data)
    df = pd.DataFrame(compound_data, columns=["!ID", "!Identifiers"])

    compound_table = SBtabTable.from_data_frame(
        df, table_id="Compound", table_type="Compound"
    )
    return compound_table,pathway_compound_all

def flux(pathway_df):
    # flux
    flux_data = []
    for ri in pathway_df.index:
        flux_tuple = ()
        flux_tuple = flux_tuple + ('rate of reaction',ri,abs(pathway_df.loc[ri,'Fluxes']))
        flux_data.append(flux_tuple)

    df = pd.DataFrame(flux_data, columns=["!QuantityType", "!Reaction","!Value"])

    flux_table = SBtabTable.from_data_frame(
        df, table_id="Flux", table_type='Quantity',unit="mM/s")
    return flux_table

def consentration(pathway_compound_all,substrateid):
    # consentration
    concentration_data = []
    for mi in pathway_compound_all:
        con_tuple = ()
        if mi == 'pi':
            con_min = 10
            con_max = 10
            con_tuple = con_tuple + ('concentration',mi,con_min,con_max)
            concentration_data.append(con_tuple)
        elif mi == 'ppi':
            con_min = 1
            con_max = 1
            con_tuple = con_tuple + ('concentration',mi,con_min,con_max)
            concentration_data.append(con_tuple)
        elif mi == 'coa':
            con_min = 1
            con_max = 1
            con_tuple = con_tuple + ('concentration',mi,con_min,con_max)
            concentration_data.append(con_tuple)
        elif mi == 'co2':
            con_min = 0.01
            con_max = 0.01
            con_tuple = con_tuple + ('concentration',mi,con_min,con_max)
            concentration_data.append(con_tuple)
        elif mi == 'atp':
            con_min = 1
            con_max = 1
            con_tuple = con_tuple + ('concentration',mi,con_min,con_max)
            concentration_data.append(con_tuple)
        elif mi == 'adp':
            con_min = 0.1
            con_max = 0.1
            con_tuple = con_tuple + ('concentration',mi,con_min,con_max)
            concentration_data.append(con_tuple)
        elif mi == 'amp':
            con_min = 0.1
            con_max = 0.1
            con_tuple = con_tuple + ('concentration',mi,con_min,con_max)
            concentration_data.append(con_tuple)
        elif mi == 'nad':
            con_min = 0.1
            con_max = 0.1
            con_tuple = con_tuple + ('concentration',mi,con_min,con_max)
            concentration_data.append(con_tuple)
        elif mi == 'nadh':
            con_min = 1
            con_max = 1
            con_tuple = con_tuple + ('concentration',mi,con_min,con_max)
            concentration_data.append(con_tuple)
        elif mi == 'nadp':
            con_min = 0.1
            con_max = 0.1
            con_tuple = con_tuple + ('concentration',mi,con_min,con_max)
            concentration_data.append(con_tuple)
        
        elif mi == 'nadph':
            con_min = 1
            con_max = 1
            con_tuple = con_tuple + ('concentration',mi,con_min,con_max)
            concentration_data.append(con_tuple)
        elif mi == 'h':
            con_min = 1
            con_max = 1
            con_tuple = con_tuple + ('concentration',mi,con_min,con_max)
            concentration_data.append(con_tuple)
        elif mi == 'h2o':
            con_min = 1
            con_max = 1
            con_tuple = con_tuple + ('concentration',mi,con_min,con_max)
            concentration_data.append(con_tuple)
        elif mi == substrateid:
            con_min = 10
            con_max = 10
            con_tuple = con_tuple + ('concentration',mi,con_min,con_max)
            concentration_data.append(con_tuple)
        else:
            con_min = 0.001
            con_max = 10
            con_tuple = con_tuple + ('concentration',mi,con_min,con_max)
            concentration_data.append(con_tuple)

    # print(consentration_data)
    df = pd.DataFrame(concentration_data, columns=["!QuantityType", "!Compound","!Min","!Max"])

    concentration_table = SBtabTable.from_data_frame(
        df, table_id="ConcentrationConstraint", table_type='Quantity',unit="mM")
    return concentration_table

def pathway_transfer_sbtab(substrate,pathwayfile,outfile,delmetcompartment):
    
    pathway_df = pd.read_table(pathwayfile)
    # /hpcfs/fhome/weif/work/BKMD/data/MetaCyc_met_mapping.csv
    metacyc_df = pd.read_csv('.data/MetaCyc_met_mapping.csv',index_col='MetaCyc')
    kegg_df = pd.read_csv('.data/kegg_mets_info.csv',index_col='id')
    ex_rxn_df = pd.read_csv('.data/bkmd_exchange_20240509.csv')
    trans_rxn_df = pd.read_csv('.data/bkmd_transport_20240509.csv')
    pathway_df.set_index('Reaction ID',inplace=True)
    
    # substrate = substrateid

    exchange_rxn = ex_rxn_df['exchange_rxn_id'].to_list()
    transport_rxn = trans_rxn_df['transport_rxn_id'].to_list()
    
    no_equi_rxn = ['ADD_h','ADD_hpc','THD2pp','ADD_respiratory1','ADD_respiratory2','ADD_respiratory3',]
    for ri in pathway_df.index: # Delete exchange and transport reactions in the pathway
        if ri in exchange_rxn or ri in transport_rxn or ri.startswith('SK') or ri.startswith('DM') or ri.startswith('R_EX') or ri in no_equi_rxn:
            pathway_df.drop(ri,inplace=True)
    pathway_df['Equation(ID)'] = pathway_df['Equation(ID)'].apply(lambda x: x.replace('-->','<=>').replace('<--', '<=>')) 
    
    config_table = config()
    reaction_table = reaction(pathway_df,metacyc_df,delmetcompartment)
    compound_table,pathway_compound_all = compound(pathway_df,metacyc_df,kegg_df,delmetcompartment)
    flux_table = flux(pathway_df)
    concentration_table = consentration(pathway_compound_all,substrate)

    sbtab_doc = SBtabDocument(name='pathway sbtab')
    sbtab_doc.add_sbtab(config_table)
    sbtab_doc.add_sbtab(reaction_table)
    sbtab_doc.add_sbtab(compound_table)
    sbtab_doc.add_sbtab(flux_table)
    sbtab_doc.add_sbtab(concentration_table)
    sbtab_doc.write(outfile)

def pathway_transfer_sbtab_folder(substrate,product,delmetcompartment,outputpath,model_id=None):

    if model_id is not None:
        folder_name = f"{substrate}-{product}-{model_id}-sbtab"
    else:
        folder_name = f"{substrate}-{product}-sbtab"
    
    outfolder_equ = os.path.join(outputpath, folder_name)
    isExists = os.path.exists(outfolder_equ)  # Determine whether the path exists
    if not isExists:
        os.makedirs(outfolder_equ)
    else:
        pass

    if model_id is not None:
        pathwayfolder = f"{substrate}-{product}-{model_id}"
    else:
        pathwayfolder = f"{substrate}-{product}"
    outfolder_pathway = os.path.join(outputpath, pathwayfolder)
    for file_name in os.listdir(outfolder_pathway):
        print(file_name)
        pathwayfile = os.path.join(outfolder_pathway,file_name)
        outfile = os.path.join(outfolder_equ,file_name.split('.')[0]+'-sbtab.tsv')
        pathway_transfer_sbtab(substrate,pathwayfile,outfile,delmetcompartment)

def cal_pathway_dg(pathway_sbtab_file,sta_figurefile,mdf_figurefile,alldgfile):
    comp_contrib = ComponentContribution()
    # pp = ThermodynamicModel.from_sbtab("result/3hpp/3hpp_bkmd_equilibrator_pathway_data_try.tsv", comp_contrib=comp_contrib)
    pp = ThermodynamicModel.from_sbtab(pathway_sbtab_file, comp_contrib=comp_contrib)
    pp.update_standard_dgs()
    pp.dg_confidence = 0.95
    mdf_sol = pp.mdf_analysis()
    # plt.style.use("seaborn-dark")
    plt.style.use("seaborn-v0_8-dark")
    
    # Modification of mdf_solution.py file
    fig, ax = plt.subplots(1, 1, figsize=(25, 7), dpi=300)
    detaG = mdf_sol.plot_standard_dg(ax=ax)
    ax.axes.xaxis.grid(True, which="major")
    plt.savefig(sta_figurefile,bbox_inches="tight")

    # Draw MDF diagram
    fig, ax = plt.subplots(1, 1, figsize=(25, 7), dpi=300)
    MDF = mdf_sol.plot_driving_forces(ax=ax)
    ax.axes.xaxis.grid(True, which="major")
    plt.savefig(mdf_figurefile,bbox_inches="tight")

    mdf_sol.reaction_df.to_csv(alldgfile, sep='\t', index=False) # Alldgfile tsv file
    return detaG,MDF

def round_with_units(value):
    try:
        number, unit = value.split(' ', 1)
        rounded_number = round(float(number), 2)
        return f"{rounded_number} {unit}"
    except Exception as e:
        print(f"Error processing value: {value}, type: {type(value)}, error: {e}")
        return value


def cal_pathway_dg_data(pathway_sbtab_file,sta_datafile,mdf_datafile,alldgfile):
    # pathway_sbtab_file = '/home/wei_f/QPATH_website/result/pathway2/glc__D_c-3hpp_c-iML1515-sbtab/3hpp_c-0-sbtab.tsv'
    comp_contrib = ComponentContribution()
    # pp = ThermodynamicModel.from_sbtab("result/3hpp/3hpp_bkmd_equilibrator_pathway_data_try.tsv", comp_contrib=comp_contrib)
    pp = ThermodynamicModel.from_sbtab(pathway_sbtab_file, comp_contrib=comp_contrib)
    pp.update_standard_dgs()
    pp.dg_confidence = 0.95
    mdf_sol = pp.mdf_analysis()

    data_df = mdf_sol.reaction_df.copy()
    data_df.reindex()

    data_df["cml_dgm"] = mdf_sol.reaction_df.physiological_dg_prime.cumsum().apply(
        lambda x: x.m_as("kJ/mol")
    )
    data_df["cml_dg_opt"] = mdf_sol.reaction_df.optimized_dg_prime.cumsum().apply(
        lambda x: x.m_as("kJ/mol")
    )
    data_df["cml_dg_sta"] = mdf_sol.reaction_df.standard_dg_prime.cumsum().apply( 
                lambda x: x.m_as("kJ/mol")
    )

    data_df
    # xticks = 0.5 + np.arange(data_df.shape[0])
    xticks = 1 + np.arange(data_df.shape[0]) 
    xticklabels = ['0'] + data_df.reaction_id.tolist()


    yvalues_phy = [0.0] + data_df.cml_dgm.tolist()
    yvalues_opt = [0.0] + data_df.cml_dg_opt.tolist()
    yvalues_sta = [0.0] + data_df.cml_dg_sta.tolist()
    yvalues_show_opt = data_df['optimized_dg_prime'].apply(lambda x: x.m_as("kJ/mol")).tolist()
    yvalues_show_opt = [round(value, 2) for value in yvalues_show_opt]

    yvalues_show_sta = data_df['standard_dg_prime'].apply(lambda x: x.m_as("kJ/mol")).tolist()
    yvalues_show_sta = [round(value, 2) for value in yvalues_show_sta]

    sta_data = dict(zip(xticklabels, yvalues_sta)) # Generate a dictionary from two lists
    phy_data = dict(zip(xticklabels, yvalues_phy))
    opt_data = dict(zip(xticklabels, yvalues_opt))
    show_data_opt = dict(zip(xticklabels[1:], yvalues_show_opt))
    show_data_sta = dict(zip(xticklabels[1:], yvalues_show_sta))

    bottleneck_data = {}
    for i in data_df[data_df.shadow_price != 0].index:
        bottle_id = data_df.loc[i,'reaction_id']
        index = xticklabels.index(bottle_id)
        bottleneck_data[xticklabels[index - 1]] = opt_data[xticklabels[index - 1]]
        bottleneck_data[xticklabels[index]] = opt_data[xticklabels[index]]

    coll = ['flux','original_standard_dg_prime','standard_dg_prime','physiological_dg_prime','optimized_dg_prime']
    alldg_df =  mdf_sol.reaction_df.copy()
    # Ensure that these columns are of string type
    alldg_df[coll] = alldg_df[coll].astype(str) 
    alldg_df[coll] = alldg_df[coll].applymap(round_with_units)
    alldg_df.to_csv(alldgfile, sep='\t', index=False) # Alldgfile tsv file
    # mdf_sol.reaction_df.to_csv(alldgfile, sep='\t', index=False) # Alldgfile tsv file

    sta_data_all = {'title':round(yvalues_sta[-1],2)}
    sta_list = []
    for k,v in sta_data.items():
        stadic = {}
        stadic['x'] = k
        stadic['y'] = round(v,2)
        stadic['type'] = 'Standard concentrations'

        if k in list(show_data_sta.keys()):
            stadic['lable'] = show_data_sta[k]
        else:
            stadic['lable'] = 'none'
        
        sta_list.append(stadic)

    sta_data_all['data'] = sta_list

    with open(sta_datafile, 'w') as json_file:
        json.dump(sta_data_all,json_file,indent=4, ensure_ascii=False)

    # MDF JSON
    mdf_data_all = {'title':round(mdf_sol.score,2)}
    mdf_list = []

    # phy
    for k,v in phy_data.items():
        mdfdic = {}
        mdfdic['x'] = k
        mdfdic['y'] = round(v,2)
        mdfdic['type'] = 'Physiological concentrations (1 mM)'
        mdf_list.append(mdfdic)

    # mdf
    for k,v in opt_data.items():
        mdfdic = {}
        mdfdic['x'] = k
        mdfdic['y'] = round(v,2)
        mdfdic['type'] = 'MDF-optimized concentrations'

        if k in list(show_data_opt.keys()):
            mdfdic['lable'] = show_data_opt[k]
        else:
            mdfdic['lable'] = 'none'
        mdf_list.append(mdfdic)

    # Bottleneck
    for k,v in bottleneck_data.items():
        mdfdic = {}
        mdfdic['x'] = k
        mdfdic['y'] = round(v,2)
        mdfdic['type'] = 'Bottleneck reactions'
        mdf_list.append(mdfdic)

    mdf_data_all['data'] = mdf_list

    with open(mdf_datafile, 'w') as json_file:
        json.dump(mdf_data_all,json_file,indent=4, ensure_ascii=False)

def cal_pathway_dg_folder(substrate,product,outputpath,model_id=None):

    if model_id is not None:
        stafoldername = f"{substrate}-{product}-{model_id}-sta-figure"
    else:
        stafoldername = f"{substrate}-{product}-sta-figure"

    outstafolder = os.path.join(outputpath, stafoldername)
    isExists = os.path.exists(outstafolder)  # Determine whether the path exists
    if not isExists:
        os.makedirs(outstafolder)
    else:
        pass
    
    if model_id is not None:
        mdffoldername = f"{substrate}-{product}-{model_id}-mdf-figure"
    else:
        mdffoldername = f"{substrate}-{product}-mdf-figure"

    outmdffolder = os.path.join(outputpath, mdffoldername)
    isExists = os.path.exists(outmdffolder)  # Determine whether the path exists
    if not isExists:
        os.makedirs(outmdffolder)
    else:
        pass
    
    if model_id is not None:
        dgfoldername = f"{substrate}-{product}-{model_id}-dg-file"
    else:
        dgfoldername = f"{substrate}-{product}-dg-file"
    outdgfolder = os.path.join(outputpath, dgfoldername)
    isExists = os.path.exists(outdgfolder)  # Determine whether the path exists
    if not isExists:
        os.makedirs(outdgfolder)
    else:
        pass

    if model_id is not None:
        inputstabfoldername = f"{substrate}-{product}-{model_id}-sbtab"
    else:
        inputstabfoldername = f"{substrate}-{product}-sbtab"

    outfolder_sbtab = os.path.join(outputpath, inputstabfoldername)
    for file_name in os.listdir(outfolder_sbtab):
        print(file_name)
        sbtabfile = os.path.join(outfolder_sbtab,file_name)
        outsta_figurefile = os.path.join(outstafolder,file_name.split('.')[0]+'-sta.png')
        outmdf_figurefile = os.path.join(outmdffolder,file_name.split('.')[0]+'-mdf.png')
        outalldgfile = os.path.join(outdgfolder,file_name.split('.')[0]+'-alldg.tsv')

        try:
            cal_pathway_dg(sbtabfile,outsta_figurefile,outmdf_figurefile,outalldgfile)
        except:
            print('error',file_name)

def cal_pathway_dg_folder_file(substrate,product,outputpath,model_id=None):

    if model_id is not None:
        stafoldername = f"{substrate}-{product}-{model_id}-sta-file"
    else:
        stafoldername = f"{substrate}-{product}-sta-file"

    outstafolder = os.path.join(outputpath, stafoldername)
    isExists = os.path.exists(outstafolder)  # Determine whether the path exists
    if not isExists:
        os.makedirs(outstafolder)
    else:
        pass
    
    if model_id is not None:
        mdffoldername = f"{substrate}-{product}-{model_id}-mdf-file"
    else:
        mdffoldername = f"{substrate}-{product}-mdf-file"

    outmdffolder = os.path.join(outputpath, mdffoldername)
    isExists = os.path.exists(outmdffolder)  # Determine whether the path exists
    if not isExists:
        os.makedirs(outmdffolder)
    else:
        pass
    
    if model_id is not None:
        dgfoldername = f"{substrate}-{product}-{model_id}-dg-file"
    else:
        dgfoldername = f"{substrate}-{product}-dg-file"
    outdgfolder = os.path.join(outputpath, dgfoldername)
    isExists = os.path.exists(outdgfolder)  # Determine whether the path exists
    if not isExists:
        os.makedirs(outdgfolder)
    else:
        pass

    if model_id is not None:
        inputstabfoldername = f"{substrate}-{product}-{model_id}-sbtab"
    else:
        inputstabfoldername = f"{substrate}-{product}-sbtab"

    outfolder_sbtab = os.path.join(outputpath, inputstabfoldername)
    for file_name in os.listdir(outfolder_sbtab):
        print(file_name)
        sbtabfile = os.path.join(outfolder_sbtab,file_name)
        outsta_figurefile = os.path.join(outstafolder,file_name.split('.')[0]+'-sta.json')
        outmdf_figurefile = os.path.join(outmdffolder,file_name.split('.')[0]+'-mdf.json')
        outalldgfile = os.path.join(outdgfolder,file_name.split('.')[0]+'-alldg.tsv')

        try:
            cal_pathway_dg_data(sbtabfile,outsta_figurefile,outmdf_figurefile,outalldgfile)
        except:
            print('error',file_name)

def bkmd_met_del_compartment(mid):
    if mid[-2] == '_':
        mid = mid[:-2]
        # met_del_com_name[met_dic['id_cal'][:-2]] = met_dic['name']
    elif mid[-3] == '_' and not mid[-1].isupper(): # Isupper() determines whether the last letter is capitalized, excluding CPD66_ _ 45 _ 45 _ 45_CCO _ 45 _ IN
        mid = mid[:-3]
    else:
        mid = mid
    return mid

