import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'root://eospublic.cern.ch//eos/opendata/cms/Run2011A/DoubleMu/AOD/12Oct2013-v1/20000/045CCED6-033F-E311-9E93-003048678F74.root'       
    )
)

process.TFileService = cms.Service ('TFileService',
    fileName = cms.string (
    'histojpsi.root'
    )
)

process.demo = cms.EDAnalyzer('JpsiAnalyzerOpen2011',
        verbose = cms.bool(True),
	# RECO Labels
		primaryVertexProducer = cms.InputTag("offlinePrimaryVertices"),
	recoMuonsLabel = cms.InputTag("muons"),   
	# RECO Configs
    	minMuPt = cms.double(2.0),# in GeV
    	maxMuEta = cms.double(2.4),
    	minMuonLeadPt = cms.double(20.0),# in GeV
	minMuonTrailPt = cms.double(4.0), # in GeV
	minJPsiMass = cms.double(2.95),# in GeV
	maxJPsiMass = cms.double(3.25)# in GeV 
)


process.p = cms.Path(process.demo)
