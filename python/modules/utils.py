import os
import sys
import math
import ROOT
import random
from numpy.polynomial import laguerre
import numpy as np


class PhysicsObject(object):
    def __init__(self, obj, pt=0., eta=0., phi=0., mass=0., keys=[]):
        self._obj = obj
        self._index = self._obj._index
        self.__dict__["pt"] = pt
        self.__dict__["eta"] = eta
        self.__dict__["phi"] = phi
        self.__dict__["mass"] = mass
        for k in keys:
            self.__dict__[k] = getattr(obj, k)

    def p4(self):
        ret = ROOT.TLorentzVector()
        ret.SetPtEtaPhiM(self.pt, self.eta, self.phi, self.mass)
        return ret

    def originalP4(self):
        return self._obj.p4()

    def __str__(self):
        return self._obj.__str__()


def deltaPhi(x, y):
    phi1 = x
    if hasattr(x,'phi'):
        phi1 = x.phi
    phi2 = y
    if hasattr(y,'phi'):
        phi2 = y.phi
    res = phi1-phi2
    while (res > math.pi):
        res -= 2 * math.pi
    while (res <= -math.pi):
        res += 2 * math.pi

    return res


def deltaR(j1, j2):
    return math.sqrt(
        (j1.eta-j2.eta)**2 +
        deltaPhi(j1.phi, j2.phi)**2
    )


def getCtauLabel(logctau):
    if logctau == 0:
        return "0"
    elif logctau > 0:
        return "1"+(("0")*logctau)
    elif logctau < 0:
        return "0p"+(("0")*(abs(logctau)-1))+"1"


def getHist(relFileName, histName):
    rootFile = ROOT.TFile(os.path.expandvars("$CMSSW_BASE/src/"+relFileName))
    hist = rootFile.Get(histName)
    if not hist:
        raise Exception("Hist file '"+histName +
                        "' not found in file '"+relFileName+"'")
    hist = hist.Clone(histName+str(random.random()))
    hist.SetDirectory(0)
    rootFile.Close()
    return hist


def getGraph(relFileName, graphName):
    rootFile = ROOT.TFile(os.path.expandvars("$CMSSW_BASE/src/"+relFileName))
    graph = rootFile.Get(graphName)
    if not graph:
        raise Exception("Graph file '"+graphName +
                        "' not found in file '"+relFileName+"'")
    graph = graph.Clone(graphName+str(random.random()))
    rootFile.Close()
    return graph


def combineHist2D(hist1, hist2, w1, w2):
    result = hist1.Clone(hist1.GetName()+hist2.GetName())
    result.SetDirectory(0)
    for ibin in range(hist1.GetNbinsX()):
        for jbin in range(hist2.GetNbinsX()):
            result.SetBinContent(ibin+1, jbin+1,
                                 w1*hist1.GetBinContent(ibin+1, jbin+1) +
                                 w2*hist2.GetBinContent(ibin+1, jbin+1)
                                 )
            result.SetBinError(ibin+1, jbin+1,
                               math.sqrt(
                                (w1*hist1.GetBinError(ibin+1, jbin+1))**2 +
                                (w2*hist2.GetBinError(ibin+1, jbin+1))**2
                                )
                               )
    return result


def getX(hist, x):
    xbin = hist.GetXaxis().FindBin(x)
    return hist.GetBinContent(xbin), hist.GetBinError(xbin)


def getSFXY(hist, x, y):
    xBin = hist.GetXaxis().FindBin(x)
    yBin = hist.GetYaxis().FindBin(y)

    if xBin == 0:
        xBin = 1
    if xBin > hist.GetNbinsX():
        xBin = hist.GetNbinsX()

    if yBin == 0:
        yBin = 1
    if yBin > hist.GetNbinsY():
        yBin = hist.GetNbinsY()

    return hist.GetBinContent(xBin, yBin), hist.GetBinError(xBin, yBin)



def getAbscissasAndWeights(N=5):
    # Laguerre polynomial roots and weights for Laguerre-Guass quadrature
    coef = np.concatenate([np.zeros(N), [1]])
    roots = laguerre.lagroots(coef)
    weights = []
    for n, root in enumerate(roots):
        n = n+1
        array = np.concatenate([np.zeros(N+1), [1]])
        value = laguerre.lagval(root, array)
        weight = root/((N+1)*value)**2
        weights.append(weight)

    return roots, weights

def calcHash(value):
    h = ((value >> 16) ^ value) * 0x45d9f3b
    h = ((h >> 16) ^ h) * 0x45d9f3b
    h = (h >> 16) ^ h
    return h
