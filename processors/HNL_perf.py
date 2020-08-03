import os
import sys
import math
import argparse
import random
import ROOT
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor \
    import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel \
    import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.modules import *

parser = argparse.ArgumentParser()

parser.add_argument('--testMode', dest='testMode', action='store_true', default=False)
parser.add_argument('--isData', dest='isData',
                    action='store_true', default=False)
parser.add_argument('--isSignal', dest='isSignal',
                    action='store_true', default=False)
parser.add_argument('--year', dest='year',
                    action='store', type=int, default=2016)
parser.add_argument('--input', dest='inputFiles', action='append', default=[])
parser.add_argument('--inv', dest='invertLeptons', action='store_true', default=False)
parser.add_argument('output', nargs=1)

args = parser.parse_args()

testMode = args.testMode
print "isData:",args.isData
print "isSignal:",args.isSignal
print "inputs:",len(args.inputFiles)

for inputFile in args.inputFiles:
    if "-2016" in inputFile or "Run2016" in inputFile:
        year = 2016
    elif "-2017" in inputFile or "Run2017" in inputFile:
        year = 2017
    elif "-2018" in inputFile or "Run2018" in inputFile:
        year = 2018
    else:
        year = args.year
    rootFile = ROOT.TFile.Open(inputFile)
    if not rootFile:
        print "CRITICAL - file '"+inputFile+"' not found!"
        sys.exit(1)
    tree = rootFile.Get("Events")
    if not tree:
        print "CRITICAL - 'Events' tree not found in file '"+inputFile+"'!"
        sys.exit(1)
    print " - ", inputFile, ", events=", tree.GetEntries()

print "year:", year
print "output directory:", args.output[0]

globalOptions = {
    "isData": args.isData,
    "year": year
}

isMC = not args.isData

minMuonPt = {2016: 25., 2017: 28., 2018: 25.}
minElectronPt = {2016: 29., 2017: 34., 2018: 34.}

if isMC:
    jecTags = {2016: 'Summer16_07Aug2017_V11_MC',
               2017: 'Fall17_17Nov2017_V32_MC',
               2018: 'Autumn18_V19_MC'
               }

    jerTags = {2016: 'Summer16_25nsV1_MC',
               2017: 'Fall17_V3_MC',
               2018: 'Autumn18_V7_MC'
               }

if args.isData:
    jecTags = {2016: 'Summer16_07Aug2017All_V11_DATA',
               2017: 'Fall17_17Nov2017_V32_DATA',
               2018: 'Autumn18_V19_DATA'
               }


leptonSelection = [
    EventSkim(selection=lambda event: event.nTrigObj > 0),
    MuonSelection(
        outputName="tightMuon",
        storeKinematics=[],
        storeWeights=True,
        muonMinPt=minMuonPt[globalOptions["year"]],
        triggerMatch=True,
        muonID=MuonSelection.TIGHT,
        muonIso=MuonSelection.INV if args.invertLeptons else MuonSelection.TIGHT,
        selectLeadingOnly=True,
        globalOptions=globalOptions
    ),
    ElectronSelection(
        outputName="tightElectron",
        storeKinematics=[],
        electronMinPt=minElectronPt[globalOptions["year"]],
        electronID="Inv" if args.invertLeptons else "Iso_WP90",
        storeWeights=True,
        triggerMatch=True,
        selectLeadingOnly=True,
        globalOptions=globalOptions
    ),
    EventSkim(selection=lambda event: (event.ntightMuon + event.ntightElectron) > 0),
    SingleMuonTriggerSelection(
        inputCollection=lambda event: event.tightMuon,
        outputName="IsoMuTrigger",
        storeWeights=True,
        globalOptions=globalOptions
    ),
    SingleElectronTriggerSelection(
        inputCollection=lambda event: event.tightElectron,
        outputName="IsoElectronTrigger",
        storeWeights=False,
        globalOptions=globalOptions
    ),
    EventSkim(selection=lambda event: (event.IsoMuTrigger_flag + event.IsoElectronTrigger_flag) > 0),
    MuonSelection(
        inputCollection=lambda event: event.tightMuon_unselected,
        outputName="looseMuons",
        storeKinematics=[],
        muonMinPt=5.,
        muonID=MuonSelection.LOOSE,
        muonIso=MuonSelection.NONE,
        globalOptions=globalOptions
    ),
    ElectronSelection(
        inputCollection=lambda event: event.tightElectron_unselected,
        outputName="looseElectrons",
        storeKinematics=[],

        electronMinPt=5.,
        electronID="Custom",
        globalOptions=globalOptions
    ),

    EventSkim(selection=lambda event: (event.ntightMuon + event.ntightElectron + event.nlooseMuons + event.nlooseElectrons ) <= 2),
    LeptonCollecting(
        tightMuonCollection=lambda event:event.tightMuon,
        tightElectronCollection=lambda event:event.tightElectron,
        looseMuonCollection=lambda event:event.looseMuons,
        looseElectronCollection=lambda event:event.looseElectrons
    )
]

analyzerChain = []

analyzerChain.extend(leptonSelection)


analyzerChain.append(
    InvariantSystem(
        inputCollection= lambda event:
            sorted(event.tightMuon+event.looseMuons+event.tightElectron+event.looseElectrons,key=lambda x: -x.pt)[:2],
        outputName="dilepton"
    )
)

featureDictFile = "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/feature_dict.py"
modelPath = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weight2016.pb",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weight2017.pb",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weight2018.pb"
}
modelGunPath = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weightGun2016.pb",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weightGun2017.pb",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weightGun2018.pb"
}

jesUncertaintyFile = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer16_07Aug2017_V11_MC_Uncertainty_AK4PFchs.txt",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.txt",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Autumn18_V19_MC_Uncertainty_AK4PFchs.txt"
}
jerResolutionFile = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Fall17_V3_MC_PtResolution_AK4PFchs.txt",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Autumn18_V7_MC_PtResolution_AK4PFchs.txt"
}

jerSFUncertaintyFile = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer16_25nsV1_MC_SF_AK4PFchs.txt",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Fall17_V3_MC_SF_AK4PFchs.txt",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Autumn18_V7_MC_SF_AK4PFchs.txt"
}


analyzerChain.append(
     MetFilter(
        globalOptions=globalOptions,
        outputName="MET_filter"
     )
)


if not isMC:
    print "ERROR - performance can only be evaluated on MC"
    sys.exit(1)


analyzerChain.append(
    JetMetUncertainties(
        jesUncertaintyFile=jesUncertaintyFile[year],
        jerResolutionFileName=jerResolutionFile[year],
        jerSFUncertaintyFileName=jerSFUncertaintyFile[year],
        jetKeys = ['pt', 'eta', 'phi' , 'jetId', 'nConstituents'],
    )
)

for systName, jetCollection in [
    ("nominal", lambda event: event.jets_nominal),
    #("jerUp", lambda event: event.jets_jerUp),
    #("jerDown", lambda event: event.jets_jerDown),
    #("jesTotalUp", lambda event: event.jets_jesTotalUp),
    #("jesTotalDown", lambda event: event.jets_jesTotalDown),
]:

    analyzerChain.append(
        JetSelection(
            inputCollection=jetCollection,
            jetMinPt=15.,
            jetId=0,
            leptonCollection=lambda event: event.leadingLepton,
            outputName="selectedJets_"+systName,
            globalOptions=globalOptions
        )
    )

    
    analyzerChain.append(
        JetTruthFlags(
            inputCollection=lambda event, systName=systName: getattr(event, "selectedJets_"+systName),
            outputName="selectedJets_"+systName,
            latentVariables=['displacement','displacement_xy','displacement_z'],
            flags={
                'isE': ['isPrompt_E'],
                'isMU': ['isPrompt_MU'],
                'isTAU': ['isPrompt_TAU'],
                
                'isB': ['isB', 'isBB', 'isGBB', 'isLeptonic_B', 'isLeptonic_C'],
                'isC': ['isC', 'isCC', 'isGCC'],
                'isUDS': ['isS', 'isUD'],
                'isG': ['isG'],
                'isPU': ['isPU'],

                'isLLP_Q': ['isLLP_Q','isLLP_QQ','isLLP_RAD'],
                
                'isLLP_E': ['isLLP_E'],
                'isLLP_MU': ['isLLP_MU'],
                'isLLP_TAU': ['isLLP_TAU'],
                
                'isLLP_QE': ['isLLP_QE','isLLP_QQE'],
                'isLLP_QMU': ['isLLP_QMU','isLLP_QQMU'],
                'isLLP_QTAU': ['isLLP_QTAU','isLLP_QQTAU'],
                
            },
            globalOptions=globalOptions
        )
        
        
    )
    
'''
analyzerChain.append(
    EventSkim(selection=lambda event: \
        getattr(event, "nselectedJets_nominal") > 0 
        #or getattr(event, "nselectedJets_jesTotalUp") > 0 \
        #or getattr(event, "nselectedJets_jesTotalDown") > 0 \
        #or getattr(event, "nselectedJets_jerUp") > 0 \
        #or getattr(event, "nselectedJets_jerDown") > 0
    )
)
'''

displacementValues = np.linspace(-3,2,5*5+1)

analyzerChain.append(
    TaggerEvaluation(
        modelPath=modelPath[year],
        featureDictFile=featureDictFile,
        inputCollections=[
            lambda event: event.selectedJets_nominal,
            #lambda event: event.selectedJets_jesTotalUp[:4],
            #lambda event: event.selectedJets_jesTotalDown[:4],
            #lambda event: event.selectedJets_jerUp[:4],
            #lambda event: event.selectedJets_jerDown[:4]
        ],
        predictionLabels=['E','MU','TAU','B','C','UDS','G','PU','LLP_Q','LLP_E','LLP_MU','LLP_TAU','LLP_QE','LLP_QMU','LLP_QTAU'],
        taggerName="llpdnnx",
        globalOptions=globalOptions,
        evalValues = displacementValues
    )
)

analyzerChain.append(
    TaggerEvaluation(
        modelPath=modelGunPath[year],
        featureDictFile=featureDictFile,
        inputCollections=[
            lambda event: event.selectedJets_nominal,
            #lambda event: event.selectedJets_jesTotalUp[:4],
            #lambda event: event.selectedJets_jesTotalDown[:4],
            #lambda event: event.selectedJets_jerUp[:4],
            #lambda event: event.selectedJets_jerDown[:4]
        ],
        predictionLabels=['E','MU','TAU','B','C','UDS','G','PU','LLP_Q','LLP_E','LLP_MU','LLP_TAU','LLP_QE','LLP_QMU','LLP_QTAU'],
        taggerName="llpdnnx_gun",
        globalOptions=globalOptions,
        evalValues = displacementValues,
    )
)


analyzerChain.append(
    JetTaggerResult(
        inputCollection=lambda event: event.selectedJets_nominal,
        taggerName="llpdnnx",
        outputName="selectedJets_nominal",
        predictionLabels=['E','MU','TAU','B','C','UDS','G','PU','LLP_Q','LLP_E','LLP_MU','LLP_TAU','LLP_QE','LLP_QMU','LLP_QTAU'],
        evalValues=displacementValues,
        formatValue=lambda num: ("%.1f"%(num)).replace(".","p").replace('-','m'),
        globalOptions={"isData": False}
    )
)

analyzerChain.append(
    JetTaggerResult(
        inputCollection=lambda event: event.selectedJets_nominal,
        taggerName="llpdnnx_gun",
        outputName="selectedJets_nominal",
        predictionLabels=['E','MU','TAU','B','C','UDS','G','PU','LLP_Q','LLP_E','LLP_MU','LLP_TAU','LLP_QE','LLP_QMU','LLP_QTAU'],
        evalValues=displacementValues,
        formatValue=lambda num: ("%.1f"%(num)).replace(".","p").replace('-','m'),
        globalOptions={"isData": False}
    )
)
'''

for systName, jetCollection in [
    ("nominal", lambda event: event.selectedJets_nominal[:4]),
    #("jerUp", lambda event: event.selectedJets_jerUp[:4]),
    #("jerDown", lambda event: event.selectedJets_jerDown[:4]),
    #("jesTotalUp", lambda event: event.selectedJets_jesTotalUp[:4]),
    #("jesTotalDown", lambda event: event.selectedJets_jesTotalDown[:4]),
]:

    analyzerChain.append(
        HNLJetSelection(
            jetCollection=jetCollection,
            promptLeptonCollection=lambda event: event.leadingLepton,
            looseLeptonsCollection=lambda event: event.subleadingLepton,
            taggerName="llpdnnx",
            jetLabels = ['LLP_Q','LLP_MU','LLP_E','LLP_TAU'],
            outputName="hnlJets_"+systName,
            globalOptions={"isData": False}
        )
    )

    
    analyzerChain.append(
       EventCategorization(
            muonsTight=lambda event: event.tightMuon, 
            electronsTight=lambda event: event.tightElectron, 
            muonsLoose=lambda event: event.looseMuons, 
            electronsLoose=lambda event: event.looseElectrons, 	
            looseLeptons=lambda event: event.subleadingLepton,
            jetsCollection=jetCollection,
            outputName="category_"+systName,
            globalOptions=globalOptions
       )
    )
'''

for systName, jetCollection, metObject in [
    ("nominal", lambda event: event.selectedJets_nominal,
        lambda event: event.met_nominal),
    #("jerUp", lambda event: event.selectedJets_jerUp,
    #    lambda event: event.met_jerUp),
    #("jerDown", lambda event: event.selectedJets_jerDown,
    #    lambda event: event.met_jerDown),
    #("jesTotalUp", lambda event: event.selectedJets_jesTotalUp,
    #    lambda event: event.met_jesTotalUp),
    #("jesTotalDown", lambda event: event.selectedJets_jesTotalDown,
    #    lambda event: event.met_jesTotalDown),
    #("unclEnUp", lambda event: event.selectedJets_nominal,
    #    lambda event: event.met_unclEnUp),
    #("unclEnDown", lambda event: event.selectedJets_nominal,
    #    lambda event: event.met_unclEnDown),
]:
    analyzerChain.extend([
        WbosonReconstruction(
            leptonCollectionName='leadingLepton',
            metObject=metObject,
            globalOptions=globalOptions,
            outputName=systName
        )
    ])

    analyzerChain.append(
        EventObservables(
            jetCollection=jetCollection,
            leptonCollection=None,
            metInput=metObject,
            globalOptions=globalOptions,
            outputName="EventObservables_"+systName
        )
    )


storeVariables = [
    [lambda tree: tree.branch("PV_npvs", "I"), lambda tree,
     event: tree.fillBranch("PV_npvs", event.PV_npvs)],
    [lambda tree: tree.branch("PV_npvsGood", "I"), lambda tree,
     event: tree.fillBranch("PV_npvsGood", event.PV_npvsGood)],
    [lambda tree: tree.branch("fixedGridRhoFastjetAll", "F"), lambda tree,
     event: tree.fillBranch("fixedGridRhoFastjetAll",
                            event.fixedGridRhoFastjetAll)],
]


if not globalOptions["isData"]:
    storeVariables.append([lambda tree: tree.branch("genweight", "F"),
                           lambda tree,
                           event: tree.fillBranch("genweight",
                           event.Generator_weight)])

    if args.isSignal:
        for coupling in range(1,68):
            storeVariables.append([
                lambda tree, coupling=coupling: tree.branch('LHEWeights_coupling_%i'%coupling,'F'),
                lambda tree, event, coupling=coupling: tree.fillBranch('LHEWeights_coupling_%i'%coupling,getattr(event,"LHEWeights_coupling_%i"%coupling)),
            ])

analyzerChain.append(EventInfo(storeVariables=storeVariables))


p = PostProcessor(
    args.output[0],
    [args.inputFiles],
    cut="((nElectron+nMuon)>0)",
    modules=analyzerChain,
    maxEvents=-1,
    friend=True
)

p.run()
