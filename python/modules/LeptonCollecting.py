import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class LeptonCollecting(Module):

    def __init__(
        self,
        tightMuonCollection = lambda event: Collection(event, "Muon"),
        tightElectronCollection = lambda event: Collection(event, "Electron"),
        looseMuonCollection = lambda event: Collection(event, "Muon"),
        looseElectronCollection = lambda event: Collection(event, "Electron"),
        outputName = "Lepton",
        globalOptions={"isData": False, "year": 2016},
        storeKinematics=["pt", "eta", "charge", "isMuon", "isElectron","relIso"]
    ):
        
        self.globalOptions = globalOptions
        self.tightMuonCollection = tightMuonCollection
        self.tightElectronCollection = tightElectronCollection
        self.looseMuonCollection = looseMuonCollection
        self.looseElectronCollection = looseElectronCollection
        self.outputName = outputName
        self.storeKinematics = storeKinematics
 
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("nleading"+self.outputName, "I")
        self.out.branch("nsubleading"+self.outputName, "I")

        for variable in self.storeKinematics:
            self.out.branch("leading"+self.outputName+"_"+variable, "F", lenVar="nleadingLepton")
            self.out.branch("subleading"+self.outputName+"_"+variable, "F", lenVar="nsubleadingLepton")


        #for variable in self.storeKinematics:
            #self.out.branch(self.outputName+"_"+variable,"F",lenVar="n"+self.outputName)
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        tightMuon = self.tightMuonCollection(event)
        tightElectron = self.tightElectronCollection(event)

        looseMuons = self.looseMuonCollection(event)
        looseElectrons = self.looseElectronCollection(event)

        for lepton in tightMuon+looseMuons:
            lepton.isMuon = 1
            lepton.isElectron = 0
            lepton.relIso = lepton.pfRelIso04_all

        for lepton in tightElectron+looseElectrons:
            lepton.isMuon = 0
            lepton.isElectron = 1
            lepton.relIso = lepton.pfRelIso03_all

        tightLepton = []
        looseLeptons = []

        tightLepton = tightMuon+tightElectron
        looseLeptons = looseMuons+looseElectrons

        tightLepton = sorted(tightLepton, key=lambda x: x.pt, reverse=True)
        # select leading only, move subleading to "loose"
        looseLeptons.extend(tightLepton[1:])
        tightLeptons = [tightLepton[0]]
        looseLeptons = sorted(looseLeptons, key=lambda x: x.pt, reverse=True)

        self.out.fillBranch("nleading"+self.outputName, len(tightLeptons))
        self.out.fillBranch("nsubleading"+self.outputName, len(looseLeptons))

        for variable in self.storeKinematics:
            self.out.fillBranch("leading"+self.outputName+"_"+variable,map(lambda lepton: getattr(lepton,variable),tightLeptons))
            self.out.fillBranch("subleading"+self.outputName+"_"+variable,map(lambda lepton: getattr(lepton,variable),looseLeptons))


        setattr(event, "leading"+self.outputName, tightLeptons)
        setattr(event, "subleading"+self.outputName, looseLeptons)

        return True
        
