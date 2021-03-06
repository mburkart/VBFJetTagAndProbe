import FWCore.ParameterSet.Config as cms

print "Running on Data"

#filter HLT paths for T&P
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt

HLTLIST = cms.VPSet(
    # VBF + Ditau trigger
    cms.PSet(
        HLT = cms.string("HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring("hltDoublePFTau20TrackPt1LooseChargedIsolationReg"),
        path2 = cms.vstring("hltMatchedVBFTwoPFJets2CrossCleanedFromDoubleLooseChargedIsoPFTau20", "hltMatchedVBFOnePFJet2CrossCleanedFromDoubleLooseChargedIsoPFTau20"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999),
    ),
)

#hltFilter = hlt.hltHighLevel.clone(
#    TriggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),
#    HLTPaths = ["HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v*", "HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v*",
#                "HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v*"],
#    andOr = cms.bool(True), # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
#    throw = cms.bool(True) #if True: throws exception if a trigger path is invalid)
#)

### ----------------------------------------------------------------------
### gen info, only from MC
### ----------------------------------------------------------------------
genInfo = cms.EDProducer("GenFiller",
        src = cms.InputTag("prunedGenParticles"),
        storeLightFlavAndGlu = cms.bool(True) # if True, store also udcs and gluons (first copy)
)

goodJets = cms.EDFilter("PATJetRefSelector",
        src = cms.InputTag("slimmedJets"),
        cut = cms.string(
                'pt > 40' # kinematics
        ),
        filter = cms.bool(True)
)

genMatchedTaus = cms.EDFilter("genMatchTauFilter",
        taus = cms.InputTag("goodTaus")
    )

goodTaus = cms.EDFilter("PATTauRefSelector",
        src = cms.InputTag("NewTauIDsEmbedded"),
        cut = cms.string(
                'abs(eta) < 2.1 ' #kinematics
                '&& abs(charge) > 0 && abs(charge) < 2 ' #sometimes 2 prongs have charge != 1
                '&& tauID("decayModeFinding") > 0.5 ' # tau ID
                '&& tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017") > 0.5 ' # tau iso - NOTE: can as well use boolean discriminators with WP
                '&& tauID("againstMuonLoose3") > 0.5 ' # anti Muon tight
                '&& tauID("againstElectronVLooseMVA6") > 0.5 ' # anti-Ele loose
        ),
        filter = cms.bool(True)
)

tightJetIdLep = cms.EDProducer("PatJetIDValueMapProducer",
        filterParams=cms.PSet(
                version = cms.string('WINTER17'),
                quality = cms.string('TIGHTLEPVETO'),
        ),
        src = cms.InputTag("slimmedJets")
)

patTriggerUnpacker = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
                                    patTriggerObjectsStandAlone = cms.InputTag("slimmedPatTrigger"),
                                    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                                    unpackFilterLabels = cms.bool(True)
)

tauPairProducer = cms.EDProducer("Tau_pair_builder",
        taus = cms.InputTag("goodTaus")
)

jetPairProducer = cms.EDProducer("Jet_pair_builder",
        jets = cms.InputTag("goodJets"),
        jetId = cms.InputTag("tightJetIdLep"),
        taus = cms.InputTag("tauPairProducer"),
        leadPtCut = cms.double(115), #actual HLT cuts
        trailPtCut = cms.double(40),
        massCut = cms.double(650)
)

overlapFilter = cms.EDFilter("Overlap_filter",
        taus = cms.InputTag("tauPairProducer"),
        jets = cms.InputTag("jetPairProducer")
)

Ntuplizer = cms.EDAnalyzer("NtuplizerVBF",
        treeName = cms.string("TagAndProbe"),
        taus = cms.InputTag("tauPairProducer"),
        jets = cms.InputTag("jetPairProducer"),
        jetId = cms.InputTag("tightJetIdLep"),
        genParticles = cms.InputTag("genInfo"),
        triggerSet = cms.InputTag("patTriggerUnpacker"),
        triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
        Vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
        puInfo = cms.InputTag("slimmedAddPileupInfo"),
        triggerListProbe = HLTLIST,
        isMC = cms.bool(False),
        filterPath = cms.string("HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg_v")
)

TAndPSeq = cms.Sequence(
#   hltFilter +
    goodJets +
    goodTaus +
#    genMatchedTaus +
    tauPairProducer +
#   tauPairFilter +
    tightJetIdLep +
    jetPairProducer +
    overlapFilter
#    genInfo
)

NtupleSeq = cms.Sequence(
    patTriggerUnpacker +
    Ntuplizer
)
