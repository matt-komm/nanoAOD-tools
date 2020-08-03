import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR, getCtauLabel


class HNLJetSelection(Module):

    def __init__(
        self,
        jetCollection=lambda event: Collection(event, "Jet"),
        promptLeptonCollection=None,
        looseLeptonsCollection=None,
        taggerName="llpdnnx",
        jetLabels = ['LLP_Q','LLP_MU','LLP_E','LLP_TAU'],
        outputName="hnlJets",
        globalOptions={"isData": False}
    ):
        self.jetCollection = jetCollection
        self.promptLeptonCollection = promptLeptonCollection
        self.looseLeptonsCollection = looseLeptonsCollection
        
        self.taggerName = taggerName
        self.jetLabels = jetLabels
        self.outputName = outputName
        self.globalOptions = globalOptions

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for label in self.jetLabels:
            self.out.branch(self.outputName+"_"+label+"_lepton_deltaR","F")
            self.out.branch(self.outputName+"_"+label+"_lepton_dxy_sig","F")
            self.out.branch(self.outputName+"_"+label+"_jet_pt","F")
            self.out.branch(self.outputName+"_"+label+"_jet_output","F")
            self.out.branch(self.outputName+"_"+label+"_jet_parameter","F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jets = self.jetCollection(event)
        
        if len(jets)==0:
            for label in self.jetLabels:
                self.out.fillBranch(self.outputName+"_"+label+"_lepton_deltaR",-1)
                self.out.fillBranch(self.outputName+"_"+label+"_lepton_dxy_sig",-1)
                self.out.fillBranch(self.outputName+"_"+label+"_lepton_deltaR",100.)
            return True
        
        promptLepton = self.promptLeptonCollection(event)
        looseLeptons = self.looseLeptonsCollection(event)
        

        bestJetsPerLabel = {}
        for label in self.jetLabels:
            bestJetsPerLabel[label] = jets[0]

        
        for jet in jets:
            taggerOutput = getattr(jet, self.taggerName)
            for label in self.jetLabels:
                if getattr(bestJetsPerLabel[label],self.taggerName)[label]['output']<taggerOutput[label]['output']:
                    bestJetsPerLabel[label] = jet

        for label in self.jetLabels:
            jet = bestJetsPerLabel[label]
            taggerResult = getattr(jet, self.taggerName)[label]


            closestLepton = None
            minDeltaR = 100.
            for looseLepton in looseLeptons:
                dR = deltaR(looseLepton,jet)
                if dR<minDeltaR:
                    minDeltaR = dR
                    closestLepton = looseLepton
                
            if closestLepton==None:
                self.out.fillBranch(self.outputName+"_"+label+"_lepton_deltaR",-1)
                self.out.fillBranch(self.outputName+"_"+label+"_lepton_dxy_sig",-1)
            else:
                self.out.fillBranch(self.outputName+"_"+label+"_lepton_deltaR",minDeltaR)
                if math.fabs(closestLepton.dxyErr)<1e-6:
                    self.out.fillBranch(self.outputName+"_"+label+"_lepton_dxy_sig",-1)
                else:
                    self.out.fillBranch(self.outputName+"_"+label+"_lepton_dxy_sig",math.fabs(closestLepton.dxy)/math.fabs(closestLepton.dxyErr))
                
                
                
            self.out.fillBranch(self.outputName+"_"+label+"_jet_pt",jet.pt)
            self.out.fillBranch(self.outputName+"_"+label+"_jet_output",taggerResult['output'])
            self.out.fillBranch(self.outputName+"_"+label+"_jet_parameter",taggerResult['parameter'])

            
        return True
