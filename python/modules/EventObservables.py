import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel \
    import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR, deltaPhi


class EventObservables(Module):

    def __init__(
        self,
        jetCollection=lambda event: Collection(event, "Jet"),
        metInput=lambda event: Object(event, "MET"),
        leptonCollection=lambda event: Collection(event, "Muon"),
        outputName="EventObservables",
        globalOptions={"isData": False}
    ):
        self.jetCollection = jetCollection
        self.metInput = metInput
        self.leptonCollection = leptonCollection
        self.outputName = outputName
        self.globalOptions = globalOptions

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch(self.outputName+"_met", "F")
        self.out.branch(self.outputName+"_met_phi", "F")
        self.out.branch(self.outputName+"_ht", "F")
        self.out.branch(self.outputName+"_mht", "F")
        self.out.branch(self.outputName+"_mht_phi", "F")
        self.out.branch(self.outputName+"_minPhi", "F")
        
        if self.leptonCollection is not None:
            self.out.branch(self.outputName+"_mT_met_Mu", "F")
            self.out.branch(self.outputName+"_met_NoMu_pt", "F")
            self.out.branch(self.outputName+"_met_NoMu_phi", "F")
            self.out.branch(self.outputName+"_mht_NoMu", "F")
            self.out.branch(self.outputName+"_mht_NoMu_phi", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        jets = self.jetCollection(event)

        met = self.metInput(event)
        self.out.fillBranch(self.outputName+"_met", met.pt)
        self.out.fillBranch(self.outputName+"_met_phi", met.phi)

        # ht and mht
        vectorSum = ROOT.TLorentzVector()
        scalarPtSum = 0.0

        for jet in jets:
            vectorSum += jet.p4()
            scalarPtSum += jet.pt

        # minPhi
        minPhi = math.pi
        for jet in jets:
            negSum = -(vectorSum-jet.p4())
            minPhi = min(minPhi, math.fabs(deltaPhi(negSum.Phi(), jet.phi)))

        self.out.fillBranch(self.outputName+"_minPhi", minPhi)
        self.out.fillBranch(self.outputName+"_ht", scalarPtSum)
        self.out.fillBranch(self.outputName+"_mht", vectorSum.Pt())
        self.out.fillBranch(self.outputName+"_mht_phi", vectorSum.Phi())
        setattr(event, self.outputName+"_minPhi", minPhi)
        setattr(event, self.outputName+"_ht", scalarPtSum)
        setattr(event, self.outputName+"_mht", vectorSum.Pt())
        setattr(event, self.outputName+"_mht_phi", vectorSum.Phi())

        if self.leptonCollection is not None:
            met_px = met.pt*math.sin(met.phi)
            met_py = met.pt*math.cos(met.phi)
            lepton = self.leptonCollection(event)
            met_px_lc = met_px + lepton.p4().Px()
            met_py_lc = met_py + lepton.p4().Py()
            met_lc = math.sqrt(met_px_lc**2 +  met_py_lc**2)
            met_l_mT = math.sqrt(2.*lepton.pt *
                                 met.pt*(1. - math.cos(lepton.phi - met.phi)))

            self.out.fillBranch(self.outputName+"_met_NoMu_pt", met_lc)
            self.out.fillBranch(self.outputName+"_mT_met_Mu", met_l_mT)
            self.out.fillBranch(self.outputName+"_mht_NoMu",
                                (vectorSum + lepton.p4()).Pt())
            setattr(event, self.outputName+"_met_NoMu_pt", met_lc)
            setattr(event, self.outputName+"_mT_met_Mu", met_l_mT)
            setattr(event, self.outputName+"_mht_NoMu",
                    (vectorSum + lepton.p4()).Pt())

        return True
