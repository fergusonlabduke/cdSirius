#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 16:18:10 2024

Functions for working with the PySirius API to Sirius 6.0.  These functions
are called by the main functions within the cdSirius node.  Compound Discoverer
results are parsed by functions from formatSpectra for input to Sirius.

@author: pleeferguson
"""

import PySirius
from PySirius import SiriusSDK
from PySirius.models.feature_import import FeatureImport
from PySirius.models.formula_candidate import FormulaCandidate
from PySirius.models.structure_candidate import StructureCandidate
import os
from io import StringIO
import time
import pandas as pd
import numpy as np
from formatSpectra import makeFeatures
from molmass import Formula



"""
Function blocks below
"""
# Expand hill notation with space between elements (CD-style)
def expandFormula(formula):
    f = Formula(formula)
    elementCounts = f.composition().dataframe()
    elementCounts.replace(1, "", inplace = True)
    elementsWithCounts = elementCounts.index+elementCounts['Count'].astype(str)
    expandedForm = elementsWithCounts.str.cat(sep = " ")
    return(expandedForm)

# Calculate mass accuracy from a formula and measured mass
def ppmAccuracy(molForm, measuredMass):
    theoreticalMass = Formula(molForm).monoisotopic_mass
    ppm = (float(measuredMass)-theoreticalMass)/theoreticalMass*1e6
    return(ppm)

# Start a local Sirius instance and connect to it
def startSirius(siriusPath, siriusUser, siriusPW):
   sdk = SiriusSDK()
   api = sdk.start_sirius(sirius_path=os.path.abspath(siriusPath), 
                          port=8080,
                          headless=True)
   #api = sdk.connect("http://localhost:8080")

   time.sleep(10)

   # Set login information
   loginCreds = {'username': siriusUser,
                 'password': siriusPW,
                 'refreshToken': None}
   loginCreds = PySirius.AccountCredentials().from_dict(loginCreds)
   # Check if user is logged-in
   isLoggedIn = api.account().is_logged_in()
   if not isLoggedIn:
       try:
           # Login into SIRIUS web services and activate default subscription if available.
           api_response = api.account().login(True, 
                                              loginCreds, 
                                              fail_when_logged_in=False, 
                                              include_subs=False)
           print(f"Sirius user {api_response.username} logged in")
       except Exception as e:
           print("Exception when calling LoginAndAccountApi->login: %s\n" % e) 
   return(api)


# Create project space
def makeProjectSpace(api, projectSpaceName, projectSpacePath):
    path = os.path.abspath(projectSpacePath)+"/"+projectSpaceName+".sirius"
    path = os.path.normpath(path)
    projectAPI = api.projects()
    ps_info = projectAPI.create_project(projectSpaceName, path_to_project=os.path.normpath(path))
    return(ps_info)
    

# Make features from a cdResult file and import into Sirius
# Creates a list of PySirius feature import objects from a 
# list of formatted feature dicts (siriusCompounds) created
# in formatSpectra.py
def importCDfeatures(cdResult, CheckedOnly, MinPeakRating, MaxMass, Limit, ps_info, api):
    siriusCompounds = makeFeatures(cdResult, 
                                   CheckedOnly, 
                                   MinPeakRating, 
                                   MaxMass, 
                                   Limit)
    feature_import =[]
    for feat in siriusCompounds:
        feature_import_from_dict = FeatureImport.from_dict(feat)
        feature_import.append(feature_import_from_dict)
    api.features().add_aligned_features(ps_info.project_id,
                                        feature_import,
                                        profile=PySirius.InstrumentProfile("ORBITRAP"), 
                                        opt_fields = ["msData"]) 


# Configure job for Sirius processing
def configureJob(profile, formulaCandidates, MS2accuracy_ppm,
                 MS1accuracy_ppm, filterByIsotopes, enforceLipidFormula,
                 performBottomUpSearch, deNovoBelowMz, formulaConstraints,
                 detectableElements, formulaSearchDBs, timeOuts,
                 doCSIFID, structureDBs, PubChemFallback, doClassyFire, 
                 doMsNovelist, msNovelistCandidates, api):
    # Load default job submission template
    jobSub = api.jobs().get_default_job_config()
    # Set parameters for job
    jobSub.spectra_search_params.enabled = False
    jobSub.formula_id_params.enabled = True
    jobSub.formula_id_params.profile = profile
    jobSub.formula_id_params.mass_accuracy_ms2ppm = float(MS2accuracy_ppm)
    jobSub.formula_id_params.filter_by_isotope_pattern = filterByIsotopes
    jobSub.formula_id_params.enforce_el_gordo_formula = enforceLipidFormula
    jobSub.formula_id_params.perform_bottom_up_search = performBottomUpSearch
    jobSub.formula_id_params.perform_denovo_below_mz = deNovoBelowMz
    jobSub.formula_id_params.enforced_formula_constraints = formulaConstraints
    jobSub.formula_id_params.detectable_elements = detectableElements
    jobSub.formula_id_params.formula_search_dbs = formulaSearchDBs
    jobSub.formula_id_params.ilp_timeout = timeOuts
    jobSub.fingerprint_prediction_params.enabled = doCSIFID
    jobSub.structure_db_search_params.enabled = doCSIFID
    jobSub.structure_db_search_params.structure_search_dbs = structureDBs
    if PubChemFallback:
        jobSub.structure_db_search_params.expansive_search_confidence_mode=api.models().ConfidenceMode.APPROXIMATE
    else:
        jobSub.structure_db_search_params.expansive_search_confidence_mode=api.models().ConfidenceMode.OFF
    jobSub.canopus_params.enabled = doClassyFire if jobSub.structure_db_search_params.enabled else False
    jobSub.ms_novelist_params.enabled = doMsNovelist if jobSub.structure_db_search_params.enabled else False
    jobSub.config_map = {'MS1MassDeviation.allowedMassDeviation': 
                         str(MS1accuracy_ppm)+' ppm'}
    return(jobSub)

# Start Sirius processing of configured job
def executeSirius(api, ps_info, jobSub):
    # Submit job
    job = api.jobs().start_job(project_id=ps_info.project_id, job_submission=jobSub)
    command = api.jobs().get_job(project_id=ps_info.project_id, job_id=job.id, opt_fields=["command"])
    print("Sirius job started with configuration: "+command.command)
    while True:
        if api.jobs().get_job(ps_info.project_id, job.id).progress.state != 'DONE':
            time.sleep(10)
        else:
            print("Sirius job completed successfully")
            break

# Get Sirius results for import to CD
def retrieveSiriusResults(api, ps_info, jobSub, saveFingerprints):
    # Get feature IDs
    featureIds = [fid.aligned_feature_id for fid in api.features().get_aligned_features(ps_info.project_id)]
    
    # Initialize lists to populate results, depending on processes performed
    features = []
    formulaResults =[] 
    if jobSub.structure_db_search_params.enabled: 
        structureResults =[]
    if jobSub.canopus_params.enabled:
        cmpdClassResults = []
    if jobSub.ms_novelist_params.enabled:
        deNovoStructureResults =[]
    
    # Retrieve results for all features 
    for featId in featureIds:
        # Get feature info with top annotation information
        feature = api.features().get_aligned_feature(ps_info.project_id, 
                                                     featId,
                                                     opt_fields=["topAnnotations"])
        
        if feature.top_annotations.formula_annotation is None:
            continue
        
        topFormulaID = feature.top_annotations.formula_annotation.formula_id
        
        # Create dict for top annotation information
        topAnnotation = {}
        # Get feature-specific information
        # Get CD compound ID value
        CDcid = feature.external_feature_id
        topAnnotation['Compounds ID'] = CDcid
        # Get feature name
        topAnnotation['Sirius Feature Name'] = feature.name
        # Parse feature neutral monoisotopic mass
        featureMass = feature.name.split("@")[0]
        # Get top formula annotation information
        topAnnotation['Sirius Top Formula'] = expandFormula(feature.top_annotations.formula_annotation.molecular_formula)
        topAnnotation['Sirius ΔMass [ppm]'] = ppmAccuracy(topAnnotation['Sirius Top Formula'], featureMass)
        topAnnotation['Sirius Top Formula Score'] = feature.top_annotations.formula_annotation.sirius_score
    
        # Get top structure annotation
        if feature.top_annotations.structure_annotation is None:
            topAnnotation['Top CSI:FingerID Name'] = ""
            topAnnotation['Top CSI:FingerID InChIKey'] = ""
            topAnnotation['Top CSI:FingerID Score'] = ""
            topAnnotation['Top CSI:FingerID Tanimoto Sim.'] = ""
            topAnnotation['Top CSI:FingerID Confid. Exact'] = ""
            topAnnotation['Top CSI:FingerID Confid. Approx.'] = ""
        else:
            topAnnotation['Top CSI:FingerID Name'] = feature.top_annotations.structure_annotation.structure_name
            topAnnotation['Top CSI:FingerID InChIKey'] = feature.top_annotations.structure_annotation.inchi_key
            topAnnotation['Top CSI:FingerID Score'] = feature.top_annotations.structure_annotation.csi_score
            topAnnotation['Top CSI:FingerID Tanimoto Sim.'] = feature.top_annotations.structure_annotation.tanimoto_similarity
            topAnnotation['Top CSI:FingerID Confid. Exact'] = feature.top_annotations.confidence_exact_match
            topAnnotation['Top CSI:FingerID Confid. Approx.'] = feature.top_annotations.confidence_approx_match
            
        # Get top ClassyFire classification
        if feature.top_annotations.compound_class_annotation is None:
            topAnnotation['Top ClassyFire Kingdom'] = ""
            topAnnotation['Top ClassyFire Superclass'] = ""
            topAnnotation['Top ClassyFire Class'] = ""
            topAnnotation['Top ClassyFire Subclass'] = ""
            topAnnotation['Top ClassyFire Level 5'] = ""
            topAnnotation['Top ClassyFire Level 6'] = ""
        else: 
            cfLineage = len(feature.top_annotations.compound_class_annotation.classy_fire_lineage)
            topAnnotation['Top ClassyFire Kingdom'] = feature.top_annotations.compound_class_annotation.classy_fire_lineage[0].name if cfLineage > 0 else ""
            topAnnotation['Top ClassyFire Superclass'] = feature.top_annotations.compound_class_annotation.classy_fire_lineage[1].name if cfLineage > 1 else ""
            topAnnotation['Top ClassyFire Class'] = feature.top_annotations.compound_class_annotation.classy_fire_lineage[2].name if cfLineage > 2 else ""
            topAnnotation['Top ClassyFire Subclass'] = feature.top_annotations.compound_class_annotation.classy_fire_lineage[3].name if cfLineage > 3 else ""
            topAnnotation['Top ClassyFire Level 5'] = feature.top_annotations.compound_class_annotation.classy_fire_lineage[4].name if cfLineage > 4 else ""
            topAnnotation['Top ClassyFire Level 6'] = feature.top_annotations.compound_class_annotation.classy_fire_lineage[5].name if cfLineage > 5 else ""

        
        # Retrieve the formula predictions with classes & fingerprints
        formulasFields = ["statistics", "compoundClasses"]
        if saveFingerprints:
            formulasFields.append("predictedFingerprint")
        formulas = api.features().get_formula_candidates(ps_info.project_id, 
                                                         featId,
                                                         opt_fields=formulasFields)
        formulas_dicts = [FormulaCandidate.to_dict(formResult) for 
                          formResult in formulas]
        
        # Get all compound class predictions if available
        for cmpdClass in formulas_dicts:
            if jobSub.canopus_params.enabled and cmpdClass['compoundClasses'] is not None:
                cmpdClass_df = pd.DataFrame.from_dict(cmpdClass['compoundClasses']['classyFireLineage'])
                cmpdClass_df.drop(columns=['type','index','levelIndex','parentId','parentName'], inplace = True)
                cmpdClass_df['siriusFormID'] = cmpdClass['formulaId']
                cmpdClass_df['Compounds ID'] = CDcid
                cmpdClass_df.rename(columns={'id': 'Class ID',
                                             'level': 'Classification Level',
                                             'name': 'Classification Name',
                                             'description': 'Description',
                                             'probability': 'Probability'}, inplace = True)
                
                # Append ClassyFire results to list
                cmpdClassResults.append(cmpdClass_df)
                
        # Get all Sirius formula predictions
        formula_df = pd.DataFrame.from_dict(formulas_dicts)
        if not formula_df.empty:
            formula_df['MS2ppmError'] = pd.DataFrame(formula_df['medianMassDeviation'].tolist())['ppm']
            
            # Assign top annotation fingerprint if requested
            if saveFingerprints:
                topFingerprint = formula_df.loc[formula_df['formulaId'] == topFormulaID, 'predictedFingerprint'].item()
                topAnnotation['topFingerprint'] = [] if topFingerprint is None else topFingerprint
            
            # Drop predictedFingerprint field from formula data frame
            formula_df.drop(['predictedFingerprint'], axis = 1, inplace = True)
            
            # Expand molecular formulas
            formula_df['Formula'] = formula_df['molecularFormula'].apply(expandFormula)
            
            # Calculate and assign ppm mass error for MS1 measurement
            formula_df['MS1 ΔMass [ppm]'] = [ppmAccuracy(x, featureMass) for x in formula_df['molecularFormula']]
            
            # Drop extraneous formula table fields
            formula_df.drop(columns=['molecularFormula',
                                     'medianMassDeviation',
                                     'compoundClasses',
                                     'siriusScoreNormalized',
                                     'zodiacScore',
                                     'fragmentationTree',
                                     'annotatedSpectrum',
                                     'isotopePatternAnnotation',
                                     'lipidAnnotation',
                                     'canopusPrediction'], inplace = True)
            
            # Assign CDcid
            formula_df['Compounds ID'] = CDcid
            
            # Rename formula table fields
            formula_df.rename(columns={'formulaId': 'siriusFormID',
                                       'adduct': 'Adduct',
                                       'rank': 'Rank',
                                       'siriusScore': 'Sirius Score',
                                       'isotopeScore': 'Isotope Score',
                                       'treeScore': 'Tree Score',
                                       'numOfExplainedPeaks': '# Explained Peaks',
                                       'numOfExplainablePeaks': '# Explainable Peaks',
                                       'totalExplainedIntensity': 'Explained Intensity',
                                       'MS2ppmError': 'Median MS2 ΔMass [ppm]'}, inplace = True)
            
            # Append Sirius formula results to list
            formulaResults.append(formula_df)
            
        
        #Append top annotation to list
        features.append(topAnnotation)
        
        # Get all CSI:FingerID database structure predictions if available
        if jobSub.structure_db_search_params.enabled:
            structures = api.features().get_structure_candidates(ps_info.project_id, 
                                                                 featId,
                                                                 opt_fields=["dbLinks"])
            structures_dicts = [StructureCandidate.to_dict(structResult) for 
                                structResult in structures]
            structure_df = pd.DataFrame.from_dict(structures_dicts)
            if not structure_df.empty:
                dbLinks = structure_df['dbLinks']
                dbLinkList = []
                for dbLink in dbLinks:
                    PubChemID = next((item for item in dbLink if item["name"] == "PUBCHEM"), None)
                    if PubChemID:
                        PubChemID = PubChemID['id']
                    DSSToxID = next((item for item in dbLink if item["name"] == "DSSTox"), None)
                    if DSSToxID:
                        DSSToxID = DSSToxID['id']
                    dbLinkList.append({"PubChem ID": PubChemID, "DSSTox ID": DSSToxID})
                dbLink_df = pd.DataFrame(dbLinkList)
                structure_df = pd.concat([structure_df, dbLink_df], axis = 1)
                
                # Calculate and assign ppm mass error for MS1 measurement
                structure_df['ΔMass [ppm]'] = [ppmAccuracy(y, featureMass) for y in structure_df['molecularFormula']]
                
                # Expand molecular formulas
                structure_df['Formula'] = structure_df['molecularFormula'].apply(expandFormula)
                
                # Drop extraneous structure table fields
                structure_df.drop(['dbLinks',
                                   'spectralLibraryMatches',
                                   'mcesDistToTopHit',
                                   'molecularFormula'], axis = 1, inplace = True)
                
                # Assign CDcid
                structure_df['Compounds ID'] = CDcid
                
                # Rename structure table fields
                structure_df.rename(columns={'formulaId': 'siriusFormID',
                                             'inchiKey': 'InChIKey',
                                             'smiles': 'SMILES',
                                             'structureName': 'Name',
                                             'xlogP': 'Log Kow',
                                             'rank': 'Rank',
                                             'csiScore': 'CSI Score',
                                             'tanimotoSimilarity': 'Tanimoto Similarity',
                                             'adduct': 'Adduct'}, inplace = True)
                
                # Append CSI:FingerID results to list
                structureResults.append(structure_df)
            
                
        # Get all MSNovelist de novo structure predictions if available
        if jobSub.ms_novelist_params.enabled:
            deNovoStructures = api.features().get_de_novo_structure_candidates(ps_info.project_id, 
                                                                 featId,
                                                                 opt_fields=["dbLinks"])
            deNovoStructures_dicts = [StructureCandidate.to_dict(deNovoResult) for 
                                      deNovoResult in deNovoStructures]
            deNovoStructure_df = pd.DataFrame.from_dict(deNovoStructures_dicts)
            if not deNovoStructure_df.empty:
                deNovodbLinks = deNovoStructure_df['dbLinks']
                deNovodbLinkList = []
                for deNovodbLink in deNovodbLinks:
                    deNovoPubChemID = next((item for item in deNovodbLink if item["name"] == "PUBCHEM"), None)
                    if deNovoPubChemID:
                        deNovoPubChemID = deNovoPubChemID['id']
                    deNovoDSSToxID = next((item for item in deNovodbLink if item["name"] == "DSSTox"), None)
                    if deNovoDSSToxID:
                        deNovoDSSToxID = deNovoDSSToxID['id']
                    deNovodbLinkList.append({"PubChem ID": deNovoPubChemID, "DSSTox ID": deNovoDSSToxID})
                deNovodbLink_df = pd.DataFrame(deNovodbLinkList)
                deNovoStructure_df = pd.concat([deNovoStructure_df, deNovodbLink_df], axis = 1)
                
                # Calculate and assign ppm mass error for MS1 measurement
                deNovoStructure_df['ΔMass [ppm]'] = [ppmAccuracy(z, featureMass) for z in deNovoStructure_df['molecularFormula']]
                
                # Expand molecular formulas
                deNovoStructure_df['Formula'] = deNovoStructure_df['molecularFormula'].apply(expandFormula)
                
                # Drop extraneous de novo structure table fields
                deNovoStructure_df.drop(['dbLinks',
                                         'spectralLibraryMatches',
                                         'molecularFormula'], axis = 1, inplace = True)
                
                # Assign CDcid
                deNovoStructure_df['Compounds ID'] = CDcid
                
                # Rename de novo structure table fields
                deNovoStructure_df.rename(columns={'formulaId': 'siriusFormID',
                                             'inchiKey': 'InChIKey',
                                             'smiles': 'SMILES',
                                             'structureName': 'Name',
                                             'xlogP': 'Log Kow',
                                             'rank': 'Rank',
                                             'csiScore': 'CSI Score',
                                             'tanimotoSimilarity': 'Tanimoto Similarity',
                                             'adduct': 'Adduct'}, inplace = True)
                
                # Append MSNovelist results to list
                deNovoStructureResults.append(deNovoStructure_df)
    
    # Assemble result tables
    results_dict = dict()
    
    # Prepare top feature annotation results
    siriusTopAnnotations = pd.DataFrame.from_dict(features)
    siriusTopAnnotations['Compounds ID'] = siriusTopAnnotations['Compounds ID'].astype('Int64')
    siriusTopAnnotations['Top CSI:FingerID Score'] = pd.to_numeric(siriusTopAnnotations['Top CSI:FingerID Score'],
                                                                   errors = 'coerce')
    siriusTopAnnotations['Top CSI:FingerID Tanimoto Sim.'] = pd.to_numeric(siriusTopAnnotations['Top CSI:FingerID Tanimoto Sim.'],
                                                                   errors = 'coerce')
    siriusTopAnnotations['Top CSI:FingerID Confid. Exact'] = pd.to_numeric(siriusTopAnnotations['Top CSI:FingerID Confid. Exact'],
                                                                   errors = 'coerce')
    siriusTopAnnotations['Top CSI:FingerID Confid. Approx.'] = pd.to_numeric(siriusTopAnnotations['Top CSI:FingerID Confid. Approx.'],
                                                                   errors = 'coerce')
    siriusTopAnnotations.replace([np.inf, -np.inf], np.nan, inplace=True)
    
    
    # Get top annotation fingerprint vectors if requested
    if saveFingerprints:
        # First get fingerprint definitions from API
        fingerprintDef = StringIO(api.projects().get_finger_id_data(ps_info.project_id, 1))
        fingerprintDef_df = pd.read_csv(fingerprintDef, sep='\t')
        # Get Sirius fingerprint predictions
        siriusFingerprints = siriusTopAnnotations['topFingerprint'].tolist()
        siriusFingerprints_df = pd.DataFrame(siriusFingerprints, index = siriusTopAnnotations['Compounds ID'])
        siriusFingerprints_df.insert(loc=0, 
                                     column= 'SiriusFeatureName',
                                     value=siriusTopAnnotations['Sirius Feature Name'].values)
        results_dict['SiriusFingerprints'] = siriusFingerprints_df
        results_dict['SiriusFingerprintDefinitions'] = fingerprintDef_df
        
        
        # Reshape top feature annotation dataframe to drop fingerprints
        siriusTopAnnotations.drop(['topFingerprint'], axis=1, inplace=True)
    
    results_dict['SiriusTopAnnotation'] = siriusTopAnnotations
    
    # Prepare Sirius formula results
    siriusFormulas = pd.concat(formulaResults, ignore_index=True)
    siriusFormulas['SiriusFormulas ID'] = siriusFormulas.index+1
    siriusFormulas['SiriusFormulas ID'] = siriusFormulas['SiriusFormulas ID'].astype('Int64')
    siriusFormulas['siriusFormID'] = siriusFormulas['siriusFormID'].astype('Int64')
    siriusFormulas['Compounds ID'] = siriusFormulas['Compounds ID'].astype('Int64')
    siriusFormulas['Rank'] = siriusFormulas['Rank'].astype('Int64')
    siriusFormulasKey = siriusFormulas[['SiriusFormulas ID', 'siriusFormID']]
    siriusFormulas.drop(['siriusFormID'], axis = 1, inplace = True)
    results_dict['SiriusFormulas'] = siriusFormulas      
    
    # Prepare CSI:FingerID results
    if jobSub.structure_db_search_params.enabled: 
        siriusStructures = pd.concat(structureResults, ignore_index=True)
        siriusStructures['SiriusStructures ID'] = siriusStructures.index+1
        siriusStructures['SiriusStructures ID'] = siriusStructures['SiriusStructures ID'].astype('Int64')
        siriusStructures['siriusFormID'] = siriusStructures['siriusFormID'].astype('Int64')
        siriusStructures['Compounds ID'] = siriusStructures['Compounds ID'].astype('Int64')
        siriusStructures['PubChem ID'] = siriusStructures['PubChem ID'].astype('Int64')
        siriusStructures['Rank'] = siriusStructures['Rank'].astype('Int64')
        siriusStructuresJoined = pd.merge(siriusStructures, siriusFormulasKey,
                                          on = 'siriusFormID',
                                          how = 'inner')
        siriusStructures = siriusStructuresJoined.drop(['siriusFormID'], axis = 1)
        results_dict['SiriusStructures'] = siriusStructures
    
    # Prepare ClassyFire results
    if jobSub.canopus_params.enabled:
        siriusCmpdClasses = pd.concat(cmpdClassResults, ignore_index = True)
        siriusCmpdClasses['SiriusClasses ID'] = siriusCmpdClasses.index+1
        siriusCmpdClasses['SiriusClasses ID'] = siriusCmpdClasses['SiriusClasses ID'].astype('Int64')
        siriusCmpdClasses['siriusFormID'] = siriusCmpdClasses['siriusFormID'].astype('Int64')
        siriusCmpdClasses['Compounds ID'] = siriusCmpdClasses['Compounds ID'].astype('Int64')
        siriusCmpdClasses['Class ID'] = siriusCmpdClasses['Class ID'].astype('Int64')
        siriusCmpdClassesJoined = pd.merge(siriusCmpdClasses, siriusFormulasKey,
                                           on = 'siriusFormID',
                                           how = 'inner')
        siriusCmpdClasses = siriusCmpdClassesJoined.drop(['siriusFormID'], axis = 1)
        results_dict['SiriusClasses'] = siriusCmpdClasses
    
    # Prepare de novo structure results
    if jobSub.ms_novelist_params.enabled:
        siriusDeNovoStructures = pd.concat(deNovoStructureResults, ignore_index=True)
        siriusDeNovoStructures['SiriusDeNovoStructures ID'] = siriusDeNovoStructures.index+1
        siriusDeNovoStructures['SiriusDeNovoStructures ID'] = siriusDeNovoStructures['SiriusDeNovoStructures ID'].astype('Int64')
        siriusDeNovoStructures['siriusFormID'] = siriusDeNovoStructures['siriusFormID'].astype('Int64')
        siriusDeNovoStructures['Compounds ID'] = siriusDeNovoStructures['Compounds ID'].astype('Int64')
        siriusDeNovoStructures['PubChem ID'] = siriusDeNovoStructures['PubChem ID'].astype('Int64')
        siriusDeNovoStructures['Rank'] = siriusDeNovoStructures['Rank'].astype('Int64')
        siriusDeNovoStructuresJoined = pd.merge(siriusDeNovoStructures, siriusFormulasKey,
                                                on = 'siriusFormID',
                                                how = 'inner')
        siriusDeNovoStructures = siriusDeNovoStructuresJoined.drop(['siriusFormID'], axis = 1)
        results_dict['SiriusDeNovoStructures'] = siriusDeNovoStructures    
    
    return(results_dict)


# Shut down API  
def shutdownSirius():
    sdk = SiriusSDK()        
    sdk.shutdown_sirius()


        

   







