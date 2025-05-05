#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 14:04:19 2024

Core functions for running Sirius within the Compound Discoverer environment.
The main function accepts arguments from the CD node output and then calls
functions within submitJob to interact with the Sirius API.  Output is parsed
to tables, which can then be imported into the CD result file.

@author: pleeferguson
"""

import sys
import os
import json
import pandas as pd
from CdScriptingNodeHelper import ScriptingResponse
from submitJob import startSirius,makeProjectSpace,importCDfeatures,configureJob,executeSirius,retrieveSiriusResults,shutdownSirius
import string
import random
import time
#from rdkit import Chem

def print_error(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    
# Write table to data file
def writeTable(outTable, outName, outPath, withIndex = False):
    outFilename = outName+'.txt'
    outFile_path = os.path.join(outPath, outFilename)
    outTable.to_csv(outFile_path, 
                    sep='\t', 
                    encoding='utf-8', 
                    index=withIndex, 
                    header=True, 
                    quoting=1,
                    na_rep ='')

def main():
    print_error('cdSirius node')
                
    # start in development mode where nodeargs are given explicitely rather than reading it as command line argument
    if sys.argv[1] == '-devel':    
        print(f'Development mode: Current Dir is {os.getcwd()}')
        nodeargs_path = 'node_args.json'        
    else:
        nodeargs_path = sys.argv[1]        

    # parse node args from Compound Discoverer
    try:
        with open(nodeargs_path, 'r') as rf:
            nodeargs = json.load(rf) 
            compounds_path = ''
            response_path = nodeargs['ExpectedResponsePath']
            cdResult_path = nodeargs['ResultFilePath']
            parameters = nodeargs['NodeParameters']
            tables = nodeargs['Tables']
            for table in tables:
                if table['TableName'] == 'Compounds':                    
                    compounds_path = table['DataFile']                                    
                    if table['DataFormat'] != 'CSV':
                        print_error(f"Unknown Data Format {table['DataFormat']}")    
                        exit(1)
            print_error('Node arguments parsed')
    except Exception as e: 
        print_error('Could not read Compound Discoverer node args')
        print_error(str(e))
        exit(1)

    if not compounds_path:
        print_error('Compounds file not defined in node args.')
        exit(1)
    
    # Parse Sirius start & login variables from parameters
    siriusPath = str(parameters['Sirius Program Path'])
    siriusUser = str(parameters['Sirius Username'])
    siriusPW = str(parameters['Sirius Password'])
    saveSirius = parameters['Save Sirius Result'] == "True"
    saveFingerprints = parameters['Save Sirius Predicted Fingerprints'] == "True"

    # Parse CD result file path, scratch folder path, generate Sirius project space name
    cdResult = cdResult_path
    if saveSirius:
        projectSpaceName = os.path.splitext(os.path.basename(cdResult_path))[0]
    else:  
        projectSpaceName = ''.join(random.choice(string.ascii_lowercase) for x in range(10))
    (workdir, _ ) = os.path.split(response_path)
    projectSpacePath = os.path.split(cdResult_path)[0] if saveSirius else workdir
    
    # Parse CD compound settings
    CheckedOnly = parameters['Checked Feature Status Handling'] == "Checked"
    MinPeakRating = float(parameters['Peak Quality Threshold'])
    MaxMass = float(parameters['Maximum MW'])
    Limit = 10000 # Arbitrary limit switch for debugging

    # Sirius job parameters    
    # Formula prediction params
    profile = "ORBITRAP"
    formulaCandidates = int(parameters['Maximum Formula Candidates'])
    MS1accuracy_ppm = float(parameters['MS1 Mass Accuracy [ppm]'])
    MS2accuracy_ppm = float(parameters['MS2 Mass Accuracy [ppm]'])
    filterByIsotopes = parameters['Filter by Isotope Pattern'] == "True"
    enforceLipidFormula = parameters['Enforce Lipid Detection Filtering'] == "True"
    performBottomUpSearch = parameters['Perform Bottom-Up Formula Search'] == "True"
    deNovoBelowMz = float(parameters['De novo Formula Generation Threshold'])
    formulaConstraints = str(parameters['Formula Elemental Constraints'])
    detectableElements = ['B', 'S', 'Cl', 'Se', 'Br']
    formulaSearchDBs = None
    timeOuts = {"numberOfSecondsPerDecomposition": 0,
                "numberOfSecondsPerInstance": 0}
    
    # CSI-FingerID params
    doCSIFID = parameters['Predict Structures']== "True"
    structureDBs = ['DSSTOX', 'PUBCHEM']
    PubChemFallback = parameters['PubChem as Fallback'] == "True"
    
    # ClassyFire params
    doClassyFire = parameters['Predict Compound Classes']== "True"
    
    # msNovelist params
    doMsNovelist = parameters['Predict de Novo Structures']== "True"
    msNovelistCandidates = int(parameters['De Novo Structure Candidates Limit'])
    
    # Check for rational Sirius settings
    if saveFingerprints and not doCSIFID:
        print_error('WARNING: CSI:FingerID must be enabled for fingerprint prediction.  Predicted fingerprints will not be saved')
        saveFingerprints = False
    if doClassyFire and not doCSIFID:
        print_error('WARNING: CSI:FingerID must be enabled for Compound Class prediction.  No Compound Class predictions will be performed')
        doClassyFire = False
    if doMsNovelist and not doCSIFID:
        print_error('WARNING: CSI:FingerID must be enabled for de Novo Structure prediction.  No de Novo structures will be predicted')
        doMsNovelist = False
                   
    # Import Compounds table from CD export   
    try:
        compoundsExport = pd.read_csv(compounds_path, 
                                      header=0,
                                      sep = '\t')

    except Exception as e:
        print_error('Could not process exported CD Compounds table')
        print_error(e)        
        exit(1)
    
    # Start the Sirius API
    try:
        api = startSirius(siriusPath, siriusUser, siriusPW)
        print_error('Sirius API started')
        
    except Exception as e:
        print_error('Could not start the Sirius engine')
        print_error(e)
        exit(1)

    # Create and connect to Sirius project space
    try:
        ps_info = makeProjectSpace(api, projectSpaceName, projectSpacePath)
        print_error(f"Sirius project space created with name {projectSpaceName}")
        
    except Exception as e:
        print_error('Could not create Sirius project space')
        print_error(e)
        shutdownSirius()
        exit(1)
    
    # Import features from Compound Discoverer to Sirius project space
    try:
        importCDfeatures(cdResult, 
                         CheckedOnly, 
                         MinPeakRating, 
                         MaxMass, 
                         Limit, 
                         ps_info, 
                         api)
        print_error("CD compounds imported to Sirius project space")
        
    except Exception as e:
        print_error('Could not import Compound Discoverer compounds to Sirius')
        print_error(e)
        shutdownSirius()
        exit(1)
        
    # Configure job for Sirius processing
    try:
        jobSub = configureJob(profile, formulaCandidates, MS2accuracy_ppm,
                         MS1accuracy_ppm, filterByIsotopes, enforceLipidFormula,
                         performBottomUpSearch, deNovoBelowMz, formulaConstraints,
                         detectableElements, formulaSearchDBs, timeOuts,
                         doCSIFID, structureDBs, PubChemFallback, doClassyFire, 
                         doMsNovelist, msNovelistCandidates, api)
        print_error("Sirius job configured")
        
    except Exception as e:
        print_error('Could not configure Sirius job input')
        print_error(e)
        shutdownSirius()
        exit(1)
        
    # Initiate Sirius job processing
    try:
        executeSirius(api, ps_info, jobSub)
        
    except Exception as e:
        print_error('Sirius processing failed')
        print_error(e)
        shutdownSirius()
        exit(1)
        
    # Import results from Sirius job completion
    try:
        results_dict = retrieveSiriusResults(api, ps_info, jobSub, saveFingerprints)
        print_error("Sirius processing results retrieved successfully")
        time.sleep(10)
         
    except Exception as e:
        print_error('Could not retrieve Sirius processing results')
        print_error(e)
        shutdownSirius()
        exit(1)
        
    # Shut down Sirius API
    shutdownSirius()
    time.sleep(10)
    
    # Save predicted fingerprints to files if requested
    if saveFingerprints:
        writeTable(results_dict['SiriusFingerprints'],
                   os.path.splitext(os.path.basename(cdResult_path))[0]+"_fingerprints",
                   os.path.split(cdResult_path)[0],
                   withIndex=True)
        writeTable(results_dict['SiriusFingerprintDefinitions'],
                   os.path.splitext(os.path.basename(cdResult_path))[0]+"_FPkey",
                   os.path.split(cdResult_path)[0])
         
    # Export Sirius result data to tables and assemble node_response JSON
    response = ScriptingResponse()
    
    # Top annotations for Compounds Table
    topAnnotations = results_dict['SiriusTopAnnotation']
    compoundsImport = compoundsExport.merge(topAnnotations,
                                     how = 'left',
                                     on = 'Compounds ID')
    writeTable(compoundsImport, "SiriusTopAnnotations", workdir)
    response.add_table('Compounds', os.path.join(workdir, 'SiriusTopAnnotations.txt'))
    response.add_column('Compounds', 'Compounds ID', 'Int', 'ID')
    response.add_column('Compounds', 'Background', 'Boolean')
    response.add_column('Compounds', 'Sirius Feature Name', 'String')
    response.set_column_option('Compounds', 'Sirius Feature Name', 'RelativePosition', '101')
    response.add_column('Compounds', 'Sirius Top Formula', 'String')
    response.set_column_option('Compounds', 'Sirius Top Formula', 'RelativePosition', '122')
    response.add_column('Compounds', 'Sirius Top Formula Score', 'Float')
    response.set_column_option('Compounds', 'Sirius Top Formula Score', 'RelativePosition', '123')
    response.set_column_option('Compounds', 'Sirius Top Formula Score', 'FormatString', 'F1')
    response.add_column('Compounds', 'Sirius ΔMass [ppm]', 'Float')
    response.set_column_option('Compounds', 'Sirius ΔMass [ppm]', 'RelativePosition', '124')
    response.set_column_option('Compounds', 'Sirius ΔMass [ppm]', 'FormatString', 'F2')
    response.add_column('Compounds', 'Top CSI:FingerID Name', 'String')
    response.set_column_option('Compounds', 'Top CSI:FingerID Name', 'RelativePosition', '401')
    response.add_column('Compounds', 'Top CSI:FingerID InChIKey', 'String')
    response.set_column_option('Compounds', 'Top CSI:FingerID InChIKey', 'RelativePosition', '402')
    response.add_column('Compounds', 'Top CSI:FingerID Score', 'Float')
    response.set_column_option('Compounds', 'Top CSI:FingerID Score', 'RelativePosition', '403')
    response.set_column_option('Compounds', 'Top CSI:FingerID Score', 'FormatString', 'F1')
    response.add_column('Compounds', 'Top CSI:FingerID Confid. Exact', 'Float')
    response.set_column_option('Compounds', 'Top CSI:FingerID Confid. Exact', 'RelativePosition', '404')
    response.set_column_option('Compounds', 'Top CSI:FingerID Confid. Exact', 'FormatString', 'F2')
    response.add_column('Compounds', 'Top CSI:FingerID Confid. Approx.', 'Float')
    response.set_column_option('Compounds', 'Top CSI:FingerID Confid. Approx.', 'RelativePosition', '405')
    response.set_column_option('Compounds', 'Top CSI:FingerID Confid. Approx.', 'FormatString', 'F2')
    response.add_column('Compounds', 'Top CSI:FingerID Tanimoto Sim.', 'Float')
    response.set_column_option('Compounds', 'Top CSI:FingerID Tanimoto Sim.', 'RelativePosition', '406')
    response.set_column_option('Compounds', 'Top CSI:FingerID Tanimoto Sim.', 'FormatString', 'F2')
    response.add_column('Compounds', 'Top ClassyFire Kingdom', 'String')
    response.set_column_option('Compounds', 'Top ClassyFire Kingdom', 'RelativePosition', '1091')
    response.add_column('Compounds', 'Top ClassyFire Superclass', 'String')
    response.set_column_option('Compounds', 'Top ClassyFire Superclass', 'RelativePosition', '1092')
    response.add_column('Compounds', 'Top ClassyFire Class', 'String')
    response.set_column_option('Compounds', 'Top ClassyFire Class', 'RelativePosition', '1093')
    response.add_column('Compounds', 'Top ClassyFire Subclass', 'String')
    response.set_column_option('Compounds', 'Top ClassyFire Subclass', 'RelativePosition', '1094')
    response.add_column('Compounds', 'Top ClassyFire Level 5', 'String')
    response.set_column_option('Compounds', 'Top ClassyFire Level 5', 'RelativePosition', '1095')
    response.add_column('Compounds', 'Top ClassyFire Level 6', 'String')
    response.set_column_option('Compounds', 'Top ClassyFire Level 6', 'RelativePosition', '1096')
    
    
    # Formula table
    formulas = results_dict['SiriusFormulas']
    writeTable(formulas.drop(['Compounds ID'], axis = 1), 
               "SiriusFormulas", workdir)
    response.add_table('SiriusFormulas', os.path.join(workdir, 'SiriusFormulas.txt'))
    response.add_column('SiriusFormulas', 'SiriusFormulas ID', 'Int', 'ID')
    response.add_column('SiriusFormulas', 'Formula', 'String')
    response.set_column_option('SiriusFormulas', 'Formula', 'RelativePosition', '10')
    response.add_column('SiriusFormulas', 'Adduct', 'String')
    response.set_column_option('SiriusFormulas', 'Adduct', 'RelativePosition', '20')
    response.add_column('SiriusFormulas', 'MS1 ΔMass [ppm]', 'Float')
    response.set_column_option('SiriusFormulas', 'MS1 ΔMass [ppm]', 'RelativePosition', '21')
    response.set_column_option('SiriusFormulas', 'MS1 ΔMass [ppm]', 'FormatString', 'F2')
    response.add_column('SiriusFormulas', 'Median MS2 ΔMass [ppm]', 'Float')
    response.set_column_option('SiriusFormulas', 'Median MS2 ΔMass [ppm]', 'RelativePosition', '22')
    response.set_column_option('SiriusFormulas', 'Median MS2 ΔMass [ppm]', 'FormatString', 'F2')
    response.add_column('SiriusFormulas', 'Rank', 'Int')
    response.set_column_option('SiriusFormulas', 'Rank', 'RelativePosition', '30')
    response.add_column('SiriusFormulas', 'Sirius Score', 'Float')
    response.set_column_option('SiriusFormulas', 'Sirius Score', 'RelativePosition', '40')
    response.set_column_option('SiriusFormulas', 'Sirius Score', 'FormatString', 'F3')
    response.add_column('SiriusFormulas', 'Isotope Score', 'Float')
    response.set_column_option('SiriusFormulas', 'Isotope Score', 'RelativePosition', '50')
    response.set_column_option('SiriusFormulas', 'Isotope Score', 'FormatString', 'F3')
    response.add_column('SiriusFormulas', 'Tree Score', 'Float')
    response.set_column_option('SiriusFormulas', 'Tree Score', 'RelativePosition', '60')
    response.set_column_option('SiriusFormulas', 'Tree Score', 'FormatString', 'F3')
    response.add_column('SiriusFormulas', '# Explained Peaks', 'Int')
    response.set_column_option('SiriusFormulas', '# Explained Peaks', 'RelativePosition', '70')
    response.add_column('SiriusFormulas', '# Explainable Peaks', 'Int')
    response.set_column_option('SiriusFormulas', '# Explainable Peaks', 'RelativePosition', '80')
    response.add_column('SiriusFormulas', 'Explained Intensity', 'Float')
    response.set_column_option('SiriusFormulas', 'Explained Intensity', 'RelativePosition', '90')
    response.set_column_option('SiriusFormulas', 'Explained Intensity', 'FormatString', 'F2') 
    # Formula to Compounds connection table
    formulas_compounds = formulas[['SiriusFormulas ID', 'Compounds ID']]
    writeTable(formulas_compounds, 'SiriusFormulas-Compounds', workdir)
    response.add_table('SiriusFormulas-Compounds', 
                       os.path.join(workdir, 'SiriusFormulas-Compounds.txt'),
                       data_format='CSVConnectionTable')
    response.set_table_option('SiriusFormulas-Compounds', 'FirstTable', 'SiriusFormulas')
    response.set_table_option('SiriusFormulas-Compounds', 'SecondTable', 'Compounds')
    response.add_column('SiriusFormulas-Compounds', 'SiriusFormulas ID', 'Int', 'ID')
    response.add_column('SiriusFormulas-Compounds', 'Compounds ID', 'Int', 'ID')
    
    if doCSIFID:
        # Structures table
        structures = results_dict['SiriusStructures']
        #structures['Structure'] = [Chem.MolToMolBlock(Chem.rdmolfiles.MolFromSmiles(m)) for 
        #                           m in structures['SMILES']]
        writeTable(structures.drop(['Compounds ID', 'SiriusFormulas ID'], axis = 1), 
                   "SiriusStructures", workdir)
        response.add_table('SiriusStructures', os.path.join(workdir, 'SiriusStructures.txt'))
        response.add_column('SiriusStructures', 'SiriusStructures ID', 'Int', 'ID')
        #response.add_column('SiriusStructures', 'Structure', 'String')
        #response.set_column_option('SiriusStructures', 'Structure', 'RelativePosition', '10')
        #response.set_column_option('SiriusStructures', 'Structure', 'SpecialDisplay', '9ACA6BD7-EB95-4F7D-A293-F18EC06D10CF')
        response.add_column('SiriusStructures', 'Name', 'String')
        response.set_column_option('SiriusStructures', 'Name', 'RelativePosition', '20')
        response.add_column('SiriusStructures', 'Formula', 'String')
        response.set_column_option('SiriusStructures', 'Formula', 'RelativePosition', '30')
        response.add_column('SiriusStructures', 'Log Kow', 'Float')
        response.set_column_option('SiriusStructures', 'Log Kow', 'RelativePosition', '31')
        response.set_column_option('SiriusStructures', 'Log Kow', 'FormatString', 'F1')
        response.add_column('SiriusStructures', 'Adduct', 'String')
        response.set_column_option('SiriusStructures', 'Adduct', 'RelativePosition', '40')
        response.add_column('SiriusStructures', 'ΔMass [ppm]', 'Float')
        response.set_column_option('SiriusStructures', 'ΔMass [ppm]', 'RelativePosition', '45')
        response.set_column_option('SiriusStructures', 'ΔMass [ppm]', 'FormatString', 'F2')
        response.add_column('SiriusStructures', 'Rank', 'Int')
        response.set_column_option('SiriusStructures', 'Rank', 'RelativePosition', '50')
        response.add_column('SiriusStructures', 'CSI Score', 'Float')
        response.set_column_option('SiriusStructures', 'CSI Score', 'RelativePosition', '60')
        response.set_column_option('SiriusStructures', 'CSI Score', 'FormatString', 'F2')
        response.add_column('SiriusStructures', 'Tanimoto Similarity', 'Float')
        response.set_column_option('SiriusStructures', 'Tanimoto Similarity', 'RelativePosition', '70')
        response.set_column_option('SiriusStructures', 'Tanimoto Similarity', 'FormatString', 'F3')
        response.add_column('SiriusStructures', 'PubChem ID', 'Int')
        response.set_column_option('SiriusStructures', 'PubChem ID', 'RelativePosition', '80')
        response.add_column('SiriusStructures', 'DSSTox ID', 'String')
        response.set_column_option('SiriusStructures', 'DSSTox ID', 'RelativePosition', '90')
        response.add_column('SiriusStructures', 'InChIKey', 'String')
        response.set_column_option('SiriusStructures', 'InChIKey', 'RelativePosition', '91')
        response.add_column('SiriusStructures', 'SMILES', 'String')
        response.set_column_option('SiriusStructures', 'SMILES', 'RelativePosition', '92')
        # Structures to Compounds connection table
        structures_compounds = structures[['SiriusStructures ID', 'Compounds ID']]
        writeTable(structures_compounds, 'SiriusStructures-Compounds', workdir)
        response.add_table('SiriusStructures-Compounds', 
                           os.path.join(workdir, 'SiriusStructures-Compounds.txt'),
                           data_format='CSVConnectionTable')
        response.set_table_option('SiriusStructures-Compounds', 'FirstTable', 'SiriusStructures')
        response.set_table_option('SiriusStructures-Compounds', 'SecondTable', 'Compounds')
        response.add_column('SiriusStructures-Compounds', 'SiriusStructures ID', 'Int', 'ID')
        response.add_column('SiriusStructures-Compounds', 'Compounds ID', 'Int', 'ID')
        # Structures to SiriusFormulas connection table
        structures_formulas = structures[['SiriusStructures ID', 'SiriusFormulas ID']]
        writeTable(structures_formulas, 'SiriusStructures-SiriusFormulas', workdir)
        response.add_table('SiriusStructures-SiriusFormulas', 
                           os.path.join(workdir, 'SiriusStructures-SiriusFormulas.txt'),
                           data_format='CSVConnectionTable')
        response.set_table_option('SiriusStructures-SiriusFormulas', 'FirstTable', 'SiriusStructures')
        response.set_table_option('SiriusStructures-SiriusFormulas', 'SecondTable', 'SiriusFormulas')
        response.add_column('SiriusStructures-SiriusFormulas', 'SiriusStructures ID', 'Int', 'ID')
        response.add_column('SiriusStructures-SiriusFormulas', 'SiriusFormulas ID', 'Int', 'ID')
                
    if doClassyFire:
        # Classes table
        classes = results_dict['SiriusClasses']
        writeTable(classes.drop(['Compounds ID', 'SiriusFormulas ID'], axis = 1), 
                   "SiriusClasses", workdir)
        response.add_table('SiriusClasses', os.path.join(workdir, 'SiriusClasses.txt'))
        response.add_column('SiriusClasses', 'SiriusClasses ID', 'Int', 'ID')
        response.add_column('SiriusClasses', 'Classification Level', 'String')
        response.set_column_option('SiriusClasses', 'Classification Level', 'RelativePosition', '10')
        response.add_column('SiriusClasses', 'Classification Name', 'String')
        response.set_column_option('SiriusClasses', 'Classification Name', 'RelativePosition', '20')
        response.add_column('SiriusClasses', 'Description', 'String')
        response.set_column_option('SiriusClasses', 'Description', 'RelativePosition', '30')
        response.add_column('SiriusClasses', 'Class ID', 'Int')
        response.set_column_option('SiriusClasses', 'Class ID', 'RelativePosition', '40')
        response.add_column('SiriusClasses', 'Probability', 'Float')
        response.set_column_option('SiriusClasses', 'Probability', 'RelativePosition', '50')
        response.set_column_option('SiriusClasses', 'Probability', 'FormatString', 'F3')
        # Classes to Compounds connection table
        classes_compounds = classes[['SiriusClasses ID', 'Compounds ID']]
        writeTable(classes_compounds, 'SiriusClasses-Compounds', workdir)
        response.add_table('SiriusClasses-Compounds', 
                           os.path.join(workdir, 'SiriusClasses-Compounds.txt'),
                           data_format='CSVConnectionTable')
        response.set_table_option('SiriusClasses-Compounds', 'FirstTable', 'SiriusClasses')
        response.set_table_option('SiriusClasses-Compounds', 'SecondTable', 'Compounds')
        response.add_column('SiriusClasses-Compounds', 'SiriusClasses ID', 'Int', 'ID')
        response.add_column('SiriusClasses-Compounds', 'Compounds ID', 'Int', 'ID')
        # Structures to SiriusFormulas connection table
        classes_formulas = classes[['SiriusClasses ID', 'SiriusFormulas ID']]
        writeTable(classes_formulas, 'SiriusClasses-SiriusFormulas', workdir)
        response.add_table('SiriusClasses-SiriusFormulas', 
                           os.path.join(workdir, 'SiriusClasses-SiriusFormulas.txt'),
                           data_format='CSVConnectionTable')
        response.set_table_option('SiriusClasses-SiriusFormulas', 'FirstTable', 'SiriusClasses')
        response.set_table_option('SiriusClasses-SiriusFormulas', 'SecondTable', 'SiriusFormulas')
        response.add_column('SiriusClasses-SiriusFormulas', 'SiriusClasses ID', 'Int', 'ID')
        response.add_column('SiriusClasses-SiriusFormulas', 'SiriusFormulas ID', 'Int', 'ID')
        
    if doMsNovelist:
        # deNovo Structures table
        deNovoStructures = results_dict['SiriusDeNovoStructures']
        #deNovoStructureResults = []
        #for m in deNovoStructures['SMILES']:
        #    try:
        #        deNovoStructureResults.append(Chem.MolToMolBlock(Chem.rdmolfiles.MolFromSmiles(m)))
        #    except Exception:
        #        deNovoStructureResults.append("")
        #deNovoStructures['Structure'] = deNovoStructureResults
        writeTable(deNovoStructures.drop(['Compounds ID', 'SiriusFormulas ID'], axis = 1), 
                   "SiriusDeNovoStructures", workdir)
        response.add_table('SiriusDeNovoStructures', os.path.join(workdir, 'SiriusDeNovoStructures.txt'))
        response.add_column('SiriusDeNovoStructures', 'SiriusDeNovoStructures ID', 'Int', 'ID')
        #response.add_column('SiriusDeNovoStructures', 'Structure', 'String')
        #response.set_column_option('SiriusDeNovoStructures', 'Structure', 'RelativePosition', '10')
        #response.set_column_option('SiriusDeNovoStructures', 'Structure', 'SpecialDisplay', '9ACA6BD7-EB95-4F7D-A293-F18EC06D10CF')
        response.add_column('SiriusDeNovoStructures', 'Name', 'String')
        response.set_column_option('SiriusDeNovoStructures', 'Name', 'RelativePosition', '20')
        response.add_column('SiriusDeNovoStructures', 'Formula', 'String')
        response.set_column_option('SiriusDeNovoStructures', 'Formula', 'RelativePosition', '30')
        response.add_column('SiriusDeNovoStructures', 'Log Kow', 'Float')
        response.set_column_option('SiriusDeNovoStructures', 'Log Kow', 'RelativePosition', '31')
        response.set_column_option('SiriusDeNovoStructures', 'Log Kow', 'FormatString', 'F1')
        response.add_column('SiriusDeNovoStructures', 'Adduct', 'String')
        response.set_column_option('SiriusDeNovoStructures', 'Adduct', 'RelativePosition', '40')
        response.add_column('SiriusDeNovoStructures', 'ΔMass [ppm]', 'Float')
        response.set_column_option('SiriusDeNovoStructures', 'MS1 ΔMass [ppm]', 'RelativePosition', '41')
        response.set_column_option('SiriusDeNovoStructures', 'MS1 ΔMass [ppm]', 'FormatString', 'F2')
        response.add_column('SiriusDeNovoStructures', 'Rank', 'Int')
        response.set_column_option('SiriusDeNovoStructures', 'Rank', 'RelativePosition', '50')
        response.add_column('SiriusDeNovoStructures', 'CSI Score', 'Float')
        response.set_column_option('SiriusDeNovoStructures', 'CSI Score', 'RelativePosition', '60')
        response.set_column_option('SiriusDeNovoStructures', 'CSI Score', 'FormatString', 'F2')
        response.add_column('SiriusDeNovoStructures', 'Tanimoto Similarity', 'Float')
        response.set_column_option('SiriusDeNovoStructures', 'Tanimoto Similarity', 'RelativePosition', '70')
        response.set_column_option('SiriusDeNovoStructures', 'Tanimoto Similarity', 'FormatString', 'F3')
        response.add_column('SiriusDeNovoStructures', 'PubChem ID', 'Int')
        response.set_column_option('SiriusDeNovoStructures', 'PubChem ID', 'RelativePosition', '80')
        response.add_column('SiriusDeNovoStructures', 'DSSTox ID', 'String')
        response.set_column_option('SiriusDeNovoStructures', 'DSSTox ID', 'RelativePosition', '90')
        response.add_column('SiriusDeNovoStructures', 'InChIKey', 'String')
        response.set_column_option('SiriusDeNovoStructures', 'InChIKey', 'RelativePosition', '91')
        response.add_column('SiriusDeNovoStructures', 'SMILES', 'String')
        response.set_column_option('SiriusDeNovoStructures', 'SMILES', 'RelativePosition', '92')
        # de Novo Structures to Compounds connection table
        deNovoStructures_compounds = deNovoStructures[['SiriusDeNovoStructures ID', 'Compounds ID']]
        writeTable(deNovoStructures_compounds, 'SiriusDeNovoStructures-Compounds', workdir)
        response.add_table('SiriusDeNovoStructures-Compounds', 
                           os.path.join(workdir, 'SiriusDeNovoStructures-Compounds.txt'),
                           data_format='CSVConnectionTable')
        response.set_table_option('SiriusDeNovoStructures-Compounds', 'FirstTable', 'SiriusDeNovoStructures')
        response.set_table_option('SiriusDeNovoStructures-Compounds', 'SecondTable', 'Compounds')
        response.add_column('SiriusDeNovoStructures-Compounds', 'SiriusDeNovoStructures ID', 'Int', 'ID')
        response.add_column('SiriusDeNovoStructures-Compounds', 'Compounds ID', 'Int', 'ID')
        # de Novo Structures to SiriusFormulas connection table
        deNovoStructures_formulas = deNovoStructures[['SiriusDeNovoStructures ID', 'SiriusFormulas ID']]
        writeTable(deNovoStructures_formulas, 'SiriusDeNovoStructures-SiriusFormulas', workdir)
        response.add_table('SiriusDeNovoStructures-SiriusFormulas', 
                           os.path.join(workdir, 'SiriusDeNovoStructures-SiriusFormulas.txt'),
                           data_format='CSVConnectionTable')
        response.set_table_option('SiriusDeNovoStructures-SiriusFormulas', 'FirstTable', 'SiriusDeNovoStructures')
        response.set_table_option('SiriusDeNovoStructures-SiriusFormulas', 'SecondTable', 'SiriusFormulas')
        response.add_column('SiriusDeNovoStructures-SiriusFormulas', 'SiriusDeNovoStructures ID', 'Int', 'ID')
        response.add_column('SiriusDeNovoStructures-SiriusFormulas', 'SiriusFormulas ID', 'Int', 'ID')
        
    
    # Save response file to disk
    response.save(response_path)




if __name__== "__main__" :
    main()
    
    
    
    
    