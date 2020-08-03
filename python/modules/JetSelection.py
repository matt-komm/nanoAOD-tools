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

    def __init__(
         self,
         inputCollection=lambda event: Collection(event, "Jet"),
         leptonCollection=lambda event: [],
         leptonFinderCollection=lambda event: [],
         outputName="selectedJets",
         jetMinPt=15.,
         jetMaxPt=-1000.,
         jetMinEta = 0.,
         jetMaxEta=2.4,
         dRCleaning=0.4,
         flagDA=False,
         storeKinematics=['pt', 'eta', 'jetId'],
         globalOptions={"isData": False},
         jetId=-1
         ):
        self.globalOptions = globalOptions
        self.inputCollection = inputCollection
        self.leptonCollection = leptonCollection
        self.leptonFinderCollection = leptonFinderCollection
        self.outputName = outputName
        self.jetMinPt = jetMinPt
        self.jetMaxPt = jetMaxPt
        self.jetMinEta = jetMinEta
        self.jetMaxEta = jetMaxEta
        self.dRCleaning = dRCleaning
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

        for variable in self.storeKinematics:
            self.out.branch(self.outputName+"_"+variable, "F", lenVar="n"+self.outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jets = self.inputCollection(event)

        selectedJets = []
        unselectedJets = []


        #print(self.outputName)

        if self.flagDA:
            flagsDA = [0.]*event.nJet

        for jet in jets:
            if jet.pt > self.jetMinPt and (self.jetMaxPt<0 or jet.pt < self.jetMaxPt) \
                and math.fabs(jet.eta) < self.jetMaxEta and (self.jetMinEta<0 or math.fabs(jet.eta) > self.jetMinEta) \
                and (jet.jetId > self.jetId):

                #note: tagger only trained for these jets
                if jet.nConstituents<4:
                    continue

                leptons = self.leptonCollection(event)
                leptonsToFind = self.leptonFinderCollection(event)

                if self.dRCleaning > 0. and leptons is not None and len(leptons) > 0:
                    mindr = min(map(lambda lepton: deltaR(lepton, jet), leptons))
                    if mindr < self.dRCleaning:
                        unselectedJets.append(jet)
                        continue

                if self.dRCleaning > 0. and leptonsToFind is not None and len(leptonsToFind) > 0:
                    mindr = min(map(lambda lepton: deltaR(lepton, jet), leptonsToFind))
                    setattr(jet, "muon_DeltaR", mindr)
                else:
                    setattr(jet, "muon_DeltaR", -1.)
                selectedJets.append(jet)
                #print(jet.pt, jet.eta)

            else:
                unselectedJets.append(jet)
                continue

            if self.flagDA:
                flagsDA[jet._index] = 1.

        if self.flagDA:
            self.out.fillBranch(self.outputName+"_forDA", flagsDA)

        #print len(selectedJets)

        for variable in self.storeKinematics:
            self.out.fillBranch(self.outputName+"_"+variable,
                                map(lambda jet: getattr(jet, variable), selectedJets))

        setattr(event, self.outputName, selectedJets)
        setattr(event, self.outputName+"_unselected", unselectedJets)

        return True
