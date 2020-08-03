import os
import sys
import math
import json
import ROOT
import random
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import getGraph, getHist, combineHist2D, getSFXY, deltaR

class ElectronSelection(Module):
    def __init__(
        self,
        inputCollection = lambda event: Collection(event, "Electron"),
        outputName = "tightElectrons",
        electronID = "Iso_WP90",
        electronMinPt = 5.,
        electronMaxEta = 2.4,
        storeKinematics=['pt','eta'],
        storeWeights=False,
        selectLeadingOnly=False,
        globalOptions={"isData":False, "year":2016},
        triggerMatch = False
    ):
        IDs = ["Iso_WP80",
               "Iso_WP90",
               "Iso_WPL",
               "noIso_WP80",
               "noIso_WP90",
               "noIso_WPL",
               "None",
               "Custom",
               "Inv"
              ]
        self.inputCollection = inputCollection
        self.outputName = outputName
        self.electronMinPt = electronMinPt
        self.electronMaxEta = electronMaxEta
        self.storeKinematics = storeKinematics
        self.storeWeights = storeWeights
        self.selectLeadingOnly = selectLeadingOnly
        self.globalOptions = globalOptions
        self.triggerMatch = triggerMatch
        if electronID not in IDs:
            print("Undefined electron ID! Choose one of the following")
            print(IDs)
            sys.exit(1)
        elif electronID == "None":
            self.storeWeights = False
            self.electronID = lambda electron: True
        elif electronID == "Custom":
            self.storeWeights = False
            self.electronID = lambda electron:  \
                                electron.GsfEleSCEtaMultiRangeCut>0 and \
                                electron.GsfEleDEtaInSeedCut>0 and \
                                electron.GsfEleDPhiInCut>0 and \
                                electron.GsfEleFull5x5SigmaIEtaIEtaCut>0 and \
                                electron.GsfEleHadronicOverEMEnergyScaledCut>0 and \
                                electron.GsfEleEInverseMinusPInverseCut>0 and \
                                electron.GsfEleConversionVetoCut>0
        elif electronID == "Inv":
            self.electronID = lambda electron: electron.mvaFall17V2noIso_WPL<1 and electron.pfRelIso03_all<0.8
            self.storeWeights = False
            
        else:
            self.electronID = lambda electron: getattr(electron, "mvaFall17V2"+electronID)

        if triggerMatch:
            self.trigger_object = lambda event: Collection(event, "TrigObj")

        id_hist_dict = {
                2016: "2016LegacyReReco_ElectronMVAREPLACE_Fall17V2.root",
                2017: "2017_ElectronMVAREPLACE.root",
                2018: "2018_ElectronMVAREPLACE.root"
        }

        id_alias_dict = {
            "Iso_WP80" : "80",
            "noIso_WP80" : "80noiso",
            "Iso_WP90" : "90",
            "noIso_WP90" : "90noiso",
        }

        #tight id efficiency
        if self.storeWeights:
            self.idHist = getHist(
                "PhysicsTools/NanoAODTools/data/electron/{}/{}".format(globalOptions["year"], id_hist_dict[globalOptions["year"]].replace("REPLACE", id_alias_dict[electronID])),
                "EGamma_SF2D"
            )

    def triggerMatched(self, electron, trigger_object):
        if self.triggerMatch:
            trig_deltaR = math.pi
            for trig_obj in trigger_object:
                if abs(trig_obj.id) != 11:
                    continue
                trig_deltaR = min(trig_deltaR, deltaR(trig_obj, electron))
            if trig_deltaR < 0.3:
                return True
            else:
                return False
        else:
            return True
 
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("n"+self.outputName, "I")
 
        if not self.globalOptions["isData"] and self.storeWeights:
            self.out.branch(self.outputName+"_weight_id_nominal","F")
            self.out.branch(self.outputName+"_weight_id_up","F")
            self.out.branch(self.outputName+"_weight_id_down","F")
            
        for variable in self.storeKinematics:
            self.out.branch(self.outputName+"_"+variable,"F",lenVar="n"+self.outputName)
            
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = self.inputCollection(event)
        muons = Collection(event, "Muon")

        if self.triggerMatch:
            trigger_object = self.trigger_object(event)
        else:
            trigger_object = None
        
        selectedElectrons = []
        unselectedElectrons = []
 
        weight_id_nominal = []
        weight_id_up = []
        weight_id_down = []
        
        for electron in electrons:
            bitmap = electron.vidNestedWPBitmap
            # decision for each cut represented by 3 bits (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
            # Electron_vidNestedWPBitmap 
            cuts = np.empty(10)

            for i in range(10):
                cuts[i] = (bitmap >> i*3) & 0x7

            electron.MinPtCut = cuts[0]
            electron.GsfEleSCEtaMultiRangeCut = cuts[1]
            electron.GsfEleDEtaInSeedCut = cuts[2]
            electron.GsfEleDPhiInCut = cuts[3]
            electron.GsfEleFull5x5SigmaIEtaIEtaCut = cuts[4]
            electron.GsfEleHadronicOverEMEnergyScaledCut = cuts[5]
            electron.GsfEleEInverseMinusPInverseCut = cuts[6]
            electron.GsfEleRelPFIsoScaledCut = cuts[7]
            electron.GsfEleConversionVetoCut = cuts[8]
            electron.GsfEleMissingHitsCut = cuts[9]

            electron.customID = electron.GsfEleSCEtaMultiRangeCut>0 and \
                                electron.GsfEleDEtaInSeedCut>0 and \
                                electron.GsfEleDPhiInCut>0 and \
                                electron.GsfEleFull5x5SigmaIEtaIEtaCut>0 and \
                                electron.GsfEleEInverseMinusPInverseCut>0 and \
                                electron.GsfEleConversionVetoCut>0
                                #electron.GsfEleHadronicOverEMEnergyScaledCut>0 and \


            if electron.pt>self.electronMinPt and math.fabs(electron.eta)<self.electronMaxEta \
                and self.electronID(electron) and self.triggerMatched(electron, trigger_object):
                if muons is not None and len(muons) > 0:
                    mindr = min(map(lambda muon: deltaR(muon, electron), muons))
                    if mindr < 0.05:
                        unselectedElectrons.append(electron)
                        continue

                selectedElectrons.append(electron)

                if self.storeWeights:
                    weight_id,weight_id_err = getSFXY(self.idHist,electron.eta,electron.pt)
                    weight_id_nominal.append(weight_id)
                    weight_id_up.append((weight_id+weight_id_err))
                    weight_id_down.append((weight_id-weight_id_err))
            else:
                unselectedElectrons.append(electron)

        if len(selectedElectrons) > 0:
            if self.selectLeadingOnly:
                unselectedElectrons.extend(selectedElectrons[1:])
                selectedElectrons = [selectedElectrons[0]]

            if not self.globalOptions["isData"] and self.storeWeights:
                weight_id_nominal = reduce(lambda x, y: x*y, weight_id_nominal)
                weight_id_up = reduce(lambda x, y: x*y, weight_id_up)
                weight_id_down = reduce(lambda x, y: x*y, weight_id_down)

        elif not self.globalOptions["isData"] and self.storeWeights:
            weight_id_nominal = 1.
            weight_id_up = 1.
            weight_id_down = 1.

        if not self.globalOptions["isData"] and self.storeWeights:
            self.out.fillBranch(self.outputName+"_weight_id_nominal", weight_id_nominal)
            self.out.fillBranch(self.outputName+"_weight_id_up", weight_id_up)
            self.out.fillBranch(self.outputName+"_weight_id_down", weight_id_down)

        self.out.fillBranch("n"+self.outputName,len(selectedElectrons))
        for variable in self.storeKinematics:
            self.out.fillBranch(self.outputName+"_"+variable,map(lambda electron: getattr(electron,variable),selectedElectrons))


        setattr(event,self.outputName,selectedElectrons)
        setattr(event,self.outputName+"_unselected",unselectedElectrons)

        return True
        
