import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR


class JetSelection(Module):

    LOOSE = 0
    TIGHT = 1
    TIGHTLEPVETO = 2
    
    def __init__(
         self,
         inputCollection=lambda event: Collection(event, "Jet"),
         leptonCollectionDRCleaning=lambda event: [],
         leptonCollectionP4Subraction=lambda event: [],
         outputName="selectedJets",
         jetMinPt=15.,
         jetMinEta =-1.,
         jetMaxEta=2.4,
         dRCleaning=0.4,
         dRP4Subtraction=0.4,
         flagDA=False,
         storeKinematics=['pt', 'eta', 'jetId'],
         globalOptions={"isData": False},
         jetId=LOOSE
     ):
        self.globalOptions = globalOptions
        self.inputCollection = inputCollection
        self.leptonCollectionDRCleaning = leptonCollectionDRCleaning
        self.leptonCollectionP4Subraction = leptonCollectionP4Subraction
        self.outputName = outputName
        self.jetMinPt = jetMinPt
        self.jetMinEta = jetMinEta
        self.jetMaxEta = jetMaxEta
        self.dRCleaning = dRCleaning
        self.dRP4Subtraction = dRP4Subtraction
        self.flagDA = flagDA
        self.storeKinematics = storeKinematics
        self.jetId = jetId

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.flagDA:
            self.out.branch(self.outputName+"_forDA", "F", lenVar="nJet")

        self.out.branch("n"+self.outputName, "I")
        for variable in self.storeKinematics:
            self.out.branch(self.outputName+"_"+variable, "F", lenVar="n"+self.outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jets = self.inputCollection(event)

        selectedJets = []
        unselectedJets = []


        leptonsForDRCleaning = self.leptonCollectionDRCleaning(event)
        leptonsForP4Subtraction = self.leptonCollectionP4Subraction(event)

        if self.flagDA:
            flagsDA = [0.]*event.nJet

        for jet in jets:
            
            if math.fabs(jet.eta) > self.jetMaxEta:
                unselectedJets.append(jet)
                continue
                
            if (self.jetMinEta>0.) and (math.fabs(jet.eta) < self.jetMinEta):
                unselectedJets.append(jet)
                continue
                
            if (jet.jetId & (1 << self.jetId)) == 0:
                unselectedJets.append(jet)
                continue
                
            #note: tagger only trained for these jets
            if jet.nConstituents<4:
                unselectedJets.append(jet)
                continue

            
            leptonP4 = ROOT.TLorentzVector(0,0,0,0)
            if self.dRP4Subtraction > 0. and len(leptonsForP4Subtraction) > 0:
                for lepton in leptonsForP4Subtraction:
                    if deltaR(lepton,jet)<self.dRP4Subtraction:
                        leptonP4 += lepton.p4()
               
            leptonPt = leptonP4.Pt()
            jetPtLeptonSubtracted = (jet.p4()-leptonP4).Pt()
            
            setattr(jet,"ptLepton",leptonPt)
            setattr(jet,"ptLeptonSubtracted",jetPtLeptonSubtracted)
            
            
            if jetPtLeptonSubtracted<self.jetMinPt:
                unselectedJets.append(jet)
                continue

                

            if self.dRCleaning > 0. and len(leptonsForDRCleaning) > 0:
                mindr = min(map(lambda lepton: deltaR(lepton, jet), leptonsForDRCleaning))
                if mindr < self.dRCleaning:
                    unselectedJets.append(jet)
                    continue


                selectedJets.append(jet)


            if self.flagDA:
                flagsDA[jet._index] = 1.

        if self.flagDA:
            self.out.fillBranch(self.outputName+"_forDA", flagsDA)

        self.out.fillBranch("n"+self.outputName, len(selectedJets))
        for variable in self.storeKinematics:
            self.out.fillBranch(self.outputName+"_"+variable,
                                map(lambda jet: getattr(jet, variable), selectedJets))

        setattr(event, self.outputName, selectedJets)
        setattr(event, self.outputName+"_unselected", unselectedJets)

        return True
