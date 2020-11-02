import os
import sys
import math
import json
import ROOT
import random
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaPhi

class WbosonReconstruction(Module):
    def __init__(self,
        leptonCollectionName = None,
        metObject =lambda event: Object(event, "MET"),
        globalOptions={"isData":False}, 
        outputName='nominal'
    ):
        self.leptonCollectionName = leptonCollectionName
        self.metObject = metObject
        self.globalOptions=globalOptions
        self.outputName=outputName
        
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self.leptonCollectionName+"_"+self.outputName+"_mtw", "F",lenVar='n'+self.leptonCollectionName)
        self.out.branch(self.leptonCollectionName+"_"+self.outputName+"_deltaPhi", "F",lenVar='n'+self.leptonCollectionName)
            
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        leptons = getattr(event,self.leptonCollectionName)
        met = self.metObject(event)
        
        mtw = np.zeros(len(leptons))
        dPhi = np.zeros(len(leptons))
        
        for i,lepton in enumerate(leptons):
            dPhi[i] = deltaPhi(lepton,met)
            mtw[i] = math.sqrt(2*lepton.pt*met.pt*(1-math.cos(dPhi[i])))
            
        self.out.fillBranch(self.leptonCollectionName+"_"+self.outputName+"_mtw", mtw)
        self.out.fillBranch(self.leptonCollectionName+"_"+self.outputName+"_deltaPhi", dPhi)
        
        
        return True
            
            
