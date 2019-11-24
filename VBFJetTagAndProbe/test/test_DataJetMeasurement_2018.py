import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import FWCore.PythonUtilities.LumiList as LumiList

process = cms.Process("TagAndProbe")

isMC = False

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

#### handling of cms line options for tier3 submission
#### the following are dummy defaults, so that one can normally use the config changing file list by hand etc.

options = VarParsing.VarParsing ('analysis')
options.register ('skipEvents',
                  -1, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Number of events to skip")
options.register ('JSONfile',
                  "", # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "JSON file (empty for no JSON)")
if isMC:
    options.outputFile = 'NTuple_MC.root'
else:
    options.outputFile = 'NTuple_Data.root'
options.inputFiles = []

#options.register('numOfThreads',
#                 10,
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.int,
#                 'Number of threads.'
#                )
#
#options.register('numOfStreams',
#                 10,
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.int,
#                 'Number of streams.'
#                )
options.parseArguments()
# START ELECTRON CUT BASED ID SECTION
#
# Set up everything that is needed to compute electron IDs and
# add the ValueMaps with ID decisions into the event data stream
#

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       era="2018-Prompt",
                       eleIDModules=[
                           "RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff",

                           "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff",
                           "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff",

                           "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff",
                           "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff",
                       ],
                       phoIDModules=[])


#START RERUNNING OF ID TRAINING
#
# set up the rerunning of the latest tau id trainings
import RecoTauTag.RecoTau.runTauIdMVA as idemb
na = idemb.TauIDEmbedder(process, cms,
        debug=True,
        updatedTauName="NewTauIDsEmbedded",
        toKeep=["2017v2", "newDM2017v2", "deepTau2017v2p1"]
)
na.runTauID()



if not isMC: # will use 80X
    from Configuration.AlCa.autoCond import autoCond
    process.GlobalTag.globaltag = '102X_dataRun2_Sep2018ABC_v2'
    # process.GlobalTag.globaltag = '102X_dataRun2_Prompt_v13'
    process.load('VBFJetTagAndProbe.VBFJetTagAndProbe.tagAndProbeJetMeasurement_2018_cff')
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
#             'file:///storage/b/akhmet/examples_files_2017/SingleMuon_ForMax.root'
#            'file:///ceph/mburkart/TestFiles/Run2017F_Tau_ReReco31Mar2018-v1.root'
#            '/store/data/Run2017F/Tau/MINIAOD/31Mar2018-v1/00000/005E78C3-5337-E811-9B1C-008CFAC91AA8.root'
            'file:///portal/ekpbms1/home/mburkart/workdir/CMSSW_9_4_8/src/TauTagAndProbe/TauTagAndProbe/test/fitter/results/pickevents.root'
        ),
    )
else:
    process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v18' #MC 25 ns miniAODv2
    process.load('VBFJetTagAndProbe.VBFJetTagAndProbe.MCanalysisJetMeasurement_cff')
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
#            'file:///storage/b/akhmet/examples_files_2017/DYJetsToLLM50_ForMax.root'
#            '/store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/005DC030-D3F4-E711-889A-02163E01A62D.root'
#            '/store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/RECOSIMstep_94X_mc2017_realistic_v10-v1/00000/705BBD66-52F2-E711-A8A1-001A648F189A.root'
             'file:///storage/9/mburkart/VBFJetTestSample/VBFJet_Pt_3200toInf_9ACFBACD-E042-E811-BB4B-008CFA1C6414.root'
#            'file:///portal/ekpbms1/home/mburkart/workdir/CMSSW_9_4_8/src/pickevents.root'
        )
    )


if options.JSONfile:
    print "Using JSON: " , options.JSONfile
    process.source.lumisToProcess = LumiList.LumiList(filename = options.JSONfile).getVLuminosityBlockRange()

if options.inputFiles:
    process.source.fileNames = cms.untracked.vstring(options.inputFiles)

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(-1)
            )

if options.maxEvents >= -1:
    process.maxEvents.input = cms.untracked.int32(options.maxEvents)
if options.skipEvents >= 0:
    process.source.skipEvents = cms.untracked.uint32(options.skipEvents)



process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
#    numberOfThreads = cms.untracked.uint32(options.numOfThreads),
#    numberOfStreams = cms.untracked.uint32(options.numOfStreams)
)

process.p = cms.Path(
    process.hltFilter +
    process.rerunMvaIsolationSequence +
    process.NewTauIDsEmbedded +
    process.goodTaus +
    process.egammaPostRecoSeq +
    process.TAndPSeq +
    process.NtupleSeq
)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# Adding ntuplizer
process.TFileService=cms.Service('TFileService',fileName=cms.string(options.outputFile))
