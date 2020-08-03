import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR, getCtauLabel


class JetTaggerResult(Module):

    def __init__(
        self,
        inputCollection=lambda event: Collection(event, "Jet"),
        taggerName="llpdnnx",
        outputName="selectedJets",
        predictionLabels=["B", "C", "UDS", "G", "PU", "isLLP_QMU_QQMU", "isLLP_Q_QQ"],
        evalValues=range(-1, 4),
        formatValue=lambda num: str(num),
        globalOptions={"isData": False}
    ):
        self.globalOptions = globalOptions
        self.taggerName = taggerName
        self.outputName = outputName
        self.inputCollection = inputCollection
        self.predictionLabels = predictionLabels
        self.evalValues = evalValues
        self.formatValue = formatValue

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for evalValue in self.evalValues:
            for label in self.predictionLabels:
                self.out.branch(self.outputName+"_"+self.taggerName+"_"+self.formatValue(evalValue)+"_"+label,"F",lenVar="n"+self.outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jets = self.inputCollection(event)

        taggerResults = {evalValue: {className: [-1.]*len(jets) for className in self.predictionLabels} for evalValue in self.evalValues}
        for ijet, jet in enumerate(jets):
            if not hasattr(jet, self.taggerName):
                print "WARNING - jet ", jet, " has no ", self.taggerName, " result stored -> skip"
                continue
            predictions = getattr(jet,self.taggerName)
            for evalValue in self.evalValues:
                for label in self.predictionLabels:
                    taggerResults[evalValue][label][ijet] = \
                        predictions[evalValue][label]
  
        for evalValue in self.evalValues:
            for label in self.predictionLabels:
                self.out.fillBranch(
                    self.outputName+"_"+self.taggerName+"_"+self.formatValue(evalValue)+"_"+label,
                    taggerResults[evalValue][label]
                )

        return True
