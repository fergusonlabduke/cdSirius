#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 13:09:10 2024

Code to extract MS and MS2 data from Compound Discoverer and put it into
a formatted dict for import into the PySirius API functionality

@author: pleeferguson
"""


import pyeds
import statistics


def makeFeatures(cdResult, CheckedOnly, MinPeakRating, MaxMass, Limit):
    with pyeds.EDS(cdResult) as eds:
       
        # define connection path and items to keep
        path = ["Compounds", "BestHitIonInstanceItem", "MassSpectrumInfoItem"]
       
        # get MS1 and MS2 from best MS2 hit and not background
        queries = {"Compounds": f"BackgroundStatus = 0 AND ExcludedBy = -1 AND MSDepth = 2 AND PeakRatingMax > {MinPeakRating}" 
                   + (' AND Checked = 1' if CheckedOnly else ''),
                   "BestHitIonInstanceItem": f"Mass < {MaxMass}",
                   "MassSpectrumInfoItem": "MassAnalyzer IN (2, 7)",}
       
        # count records
        #counts = eds.Count("Compounds", query=queries["Compounds"])
       
        # read compound data from DB
        compounds = eds.ReadHierarchy(
            path,
            queries = queries,
            orders = {"Compounds": "MaxArea"},
            descs = {"Compounds": True},
            limits = {"Compounds": Limit})
       
        # create list for Sirius input data
        siriusCompounds = []   
                
        # extract spectra and information for Sirius input for each compound
        for cmpd in compounds: 
            # read spectra from DB
            ids = [sp.IDs for hit in cmpd.Children for sp in hit.Children]
            if len(ids) == 0:
                continue
            spectra = {sp.IDs: sp for sp in eds.ReadMany("MassSpectrumItem", ids)}
            
            # Get peak width information for compound
            peakWidths = list(eds.ReadConnected("UnknownCompoundInstanceItem", parent = cmpd, properties = ['FWHM',]))
            FWHM = [peak.FWHM for peak in peakWidths]
            meanFWHMsec = statistics.fmean(FWHM)*60
            
            # get compound-specific data
            # extract ion definition from best hit MS2 information
            ionization = [bh.IonDescription for bh in cmpd.Children 
                          if bh.BestHitType.Value == 2][0]
            # check for invalid ion definitions and format adduct string
            invalidAdducts = ["2M", "+2", "+3","-2","-3", "MeOH", "ACN", "-e", "+e"]
            if any(x in ionization for x in invalidAdducts):
                if cmpd.Polarity.Value == 1:   
                    ionization = ["[M+H]+"]
                    ionMass  = cmpd.MolecularWeight + 1.00727663
                else:
                    ionization = ["[M-H']-"]
                    ionMass = cmpd.MolecularWeight - 1.00727663
            else:
                ionization = [ionization[:-1]]
                ionMass = [bh.Mass for bh in cmpd.Children 
                           if bh.BestHitType.Value == 2][0]
            
            # assemble Sirius input format dict
            siriusInput ={}
            #print(f"\n {cmpd.ID} {cmpd.MolecularWeight:.5f}@{cmpd.RetentionTime:.3f} {cmpd.Formula} {cmpd.Name}")
            siriusInput['name'] = f"{cmpd.MolecularWeight:.5f}@{cmpd.RetentionTime:.2f}"
            siriusInput['externalFeatureId'] = str(cmpd.ID)
            siriusInput['ionMass'] = ionMass
            siriusInput['charge'] = [bh.Charge for bh in cmpd.Children 
                                     if bh.BestHitType.Value == 2][0]
            siriusInput['detectedAdducts'] = ionization
            siriusInput['rtStartSeconds'] = cmpd.RetentionTime*60 - 0.5*meanFWHMsec
            siriusInput['rtEndSeconds'] = cmpd.RetentionTime*60 + 0.5*meanFWHMsec
            
            
            # select MS1 spectrum and append to dict
            MS1dict ={}
            MS1ID = [sp.IDs for hit in cmpd.Children 
                     if hit.BestHitType.Value == 1 for sp in hit.Children][0]
            MS1 = spectra[MS1ID]
            MS1spectrum = [(p.MZ, p.Intensity) for 
                   p in MS1.Spectrum.Centroids]
            # find indices of isotopic envelope range
            isoBegin = MS1spectrum.index(min(MS1spectrum, key=lambda s: 
                                             abs(s[0]-(siriusInput['ionMass']-1))))
            isoEnd = MS1spectrum.index(min(MS1spectrum, key=lambda s: 
                                           abs(s[0]-(siriusInput['ionMass']+5))))
            MS1spectrum = MS1spectrum[isoBegin:isoEnd]
            MS1dict['name'] = "MS1"
            MS1dict['msLevel'] = 1
            MS1dict['scanNumber'] = MS1.Spectrum.Header.ScanNumber
            MS1dict['peaks'] = [{'mz':p[0], 'intensity': p[1]} for 
                                p in MS1spectrum]
            siriusInput['mergedMs1'] = MS1dict
            
             
            # select best MS2 spectra and add to dict
            MS2IDs = [sp.IDs for hit in cmpd.Children for sp in hit.Children 
                      if sp.MSOrder.Value == 2]
            # extract mass spectrum
            MS2spectra = []
            for specID in MS2IDs:
                MS2dict ={}
                MS2 = spectra[specID]
                MS2dict['name'] = "MS2_"+str(MS2.Spectrum.Header.ScanNumber)
                MS2dict['msLevel'] = MS2.MSOrder.Value
                MS2dict['collisionEnergy'] = str(MS2.Spectrum.Event.ActivationEnergies[0])
                MS2dict['precursorMz'] = siriusInput['ionMass']
                MS2dict['scanNumber'] = MS2.Spectrum.Header.ScanNumber
                MS2spectrum = [(p.MZ, p.Intensity) for 
                               p in MS2.Spectrum.Centroids]
                MS2dict['peaks'] = [{'mz':p[0], 'intensity': p[1]} for 
                                    p in MS2spectrum]
                #indent = "\t"*MS2.MSOrder.Value
                #print(f"\t{indent}{spectrum.Spectrum}")
                MS2spectra.append(MS2dict)
            siriusInput['ms2Spectra'] = MS2spectra
            siriusCompounds.append(siriusInput)
        return(siriusCompounds)







    
    
    
