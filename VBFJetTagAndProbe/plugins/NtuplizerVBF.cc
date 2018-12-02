// -*- C++ -*-
//
// Package:    QCDTagAndProbe/QCDTagAndProbe
// Class:      NtuplizerVBF
// 
/**\class NtuplizerVBF NtuplizerVBF.cc QCDTagAndProbe/QCDTagAndProbe/plugins/NtuplizerVBF.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Maximilian Burkart
//         Created:  Tue, 03 Jul 2018 12:50:23 GMT
//
//

#ifndef NTUPLIZERVBF_H
#define NTUPLIZERVBF_H
// system include files
#include <memory>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>


// user include files
#include <TNtuple.h>
#include <TString.h>

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>

#include <FWCore/Framework/interface/Event.h>

#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <FWCore/ServiceRegistry/interface/Service.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <FWCore/Common/interface/TriggerNames.h> 
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
#include <HLTrigger/HLTcore/interface/HLTConfigProvider.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/Common/interface/TriggerResults.h>
#include <DataFormats/PatCandidates/interface/GenericParticle.h>

#include "tParameterSet.h"

#include <CommonTools/UtilAlgos/interface/TFileService.h>



//Set this variable to decide the number of triggers that you want to check simultaneously
#define NUMBER_OF_MAXIMUM_TRIGGERS 64
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class NtuplizerVBF : public edm::EDAnalyzer {
   public:
      explicit NtuplizerVBF(const edm::ParameterSet&);
      virtual ~NtuplizerVBF();

   private:
      virtual void beginJob();
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob();
      virtual void endRun(edm::Run const&, edm::EventSetup const&);

      void Initialize();
      bool hasFilters(const pat::TriggerObjectStandAlone&, const std::vector<std::string>&);
      int GenIndex(const pat::TauRef&, const edm::View<pat::GenericParticle>*);
      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::TauRefVector> tauTag_;
      edm::EDGetTokenT<pat::JetRefVector> leadJet_;
      edm::EDGetTokenT<edm::ValueMap<int>> jetIdTag_;
      edm::EDGetTokenT<edm::View<pat::GenericParticle>> genParticlesTag_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<std::vector<reco::Vertex>> VtxTag_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puTag_;
      edm::InputTag processName_;
      std::string treeName_;
      bool useGenMatch_;
      bool useHLTMatch_;
      bool isMC_;
      std::string filterPath_;

      TTree* tree_;
      TTree* triggerNamesTree_;
      
      // ---------------------------------------
      // variables to be filled in output tree
      ULong64_t indexevents_;
      Int_t runNumber_;
      Int_t lumi_;
      Int_t PS_column_;

      float MC_weight_;

      unsigned long tauTriggerBits_1_;
      unsigned long tauTriggerBits_woL3_1_;
      unsigned long tauTriggerBits_2_;
      unsigned long tauTriggerBits_woL3_2_;
      float tauPt_1_;
      float tauEta_1_;
      float tauPhi_1_;
      int tauCharge_1_;
      int tauDM_1_;
      int tau_genindex_1_;
      float tauTrkPt_1_;

      float tauPt_2_;
      float tauEta_2_;
      float tauPhi_2_;
      int tauCharge_2_;
      int tauDM_2_;
      int tau_genindex_2_;
      float tauTrkPt_2_;

      float leadJetPt_;
      float leadJetEta_;
      float leadJetPhi_;
      float trailJetPt_;
      float trailJetEta_;
      float trailJetPhi_;
      float thirdJetPt_;
      float thirdJetEta_;
      float thirdJetPhi_;
      bool jetPairL1Matched_;
      int nJets_;
      float Mjj_;

      bool tauByVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1_;
      bool tauByVLooseIsolationMVArun2017v2DBoldDMwLT2017_1_;
      bool tauByLooseIsolationMVArun2017v2DBoldDMwLT2017_1_;
      bool tauByMediumIsolationMVArun2017v2DBoldDMwLT2017_1_;
      bool tauByTightIsolationMVArun2017v2DBoldDMwLT2017_1_;
      bool tauByVTightIsolationMVArun2017v2DBoldDMwLT2017_1_;
      bool tauByVVTightIsolationMVArun2017v2DBoldDMwLT2017_1_;
      bool tauByVVLooseIsolationMVArun2017v2DBnewDMwLT2017_1_;
      bool tauByVLooseIsolationMVArun2017v2DBnewDMwLT2017_1_;
      bool tauByLooseIsolationMVArun2017v2DBnewDMwLT2017_1_;
      bool tauByMediumIsolationMVArun2017v2DBnewDMwLT2017_1_;
      bool tauByTightIsolationMVArun2017v2DBnewDMwLT2017_1_;
      bool tauByVTightIsolationMVArun2017v2DBnewDMwLT2017_1_;
      bool tauByVVTightIsolationMVArun2017v2DBnewDMwLT2017_1_;

      bool tauAgainstMuonLoose3_1_;
      bool tauAgainstMuonTight3_1_;
      bool tauAgainstElectronVLooseMVA6_1_;
      bool tauAgainstElectronLooseMVA6_1_;
      bool tauAgainstElectronMediumMVA6_1_;
      bool tauAgainstElectronTightMVA6_1_;
      bool tauAgainstElectronVTightMVA6_1_;

      bool tauByVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2_;
      bool tauByVLooseIsolationMVArun2017v2DBoldDMwLT2017_2_;
      bool tauByLooseIsolationMVArun2017v2DBoldDMwLT2017_2_;
      bool tauByMediumIsolationMVArun2017v2DBoldDMwLT2017_2_;
      bool tauByTightIsolationMVArun2017v2DBoldDMwLT2017_2_;
      bool tauByVTightIsolationMVArun2017v2DBoldDMwLT2017_2_;
      bool tauByVVTightIsolationMVArun2017v2DBoldDMwLT2017_2_;
      bool tauByVVLooseIsolationMVArun2017v2DBnewDMwLT2017_2_;
      bool tauByVLooseIsolationMVArun2017v2DBnewDMwLT2017_2_;
      bool tauByLooseIsolationMVArun2017v2DBnewDMwLT2017_2_;
      bool tauByMediumIsolationMVArun2017v2DBnewDMwLT2017_2_;
      bool tauByTightIsolationMVArun2017v2DBnewDMwLT2017_2_;
      bool tauByVTightIsolationMVArun2017v2DBnewDMwLT2017_2_;
      bool tauByVVTightIsolationMVArun2017v2DBnewDMwLT2017_2_;

      bool tauAgainstMuonLoose3_2_;
      bool tauAgainstMuonTight3_2_;
      bool tauAgainstElectronVLooseMVA6_2_;
      bool tauAgainstElectronLooseMVA6_2_;
      bool tauAgainstElectronMediumMVA6_2_;
      bool tauAgainstElectronTightMVA6_2_;
      bool tauAgainstElectronVTightMVA6_2_;

      int Nvtx_;
      float nTruePU_;
      UInt_t lastFilter_;



      //!Contains the parameters
      tVParameterSet parameters_;
      
      //! Maximum
      std::bitset<NUMBER_OF_MAXIMUM_TRIGGERS> tauTriggerBitSet_1_;
      std::bitset<NUMBER_OF_MAXIMUM_TRIGGERS> tauTriggerBitSet_woL3_1_;
      std::bitset<NUMBER_OF_MAXIMUM_TRIGGERS> tauTriggerBitSet_2_;
      std::bitset<NUMBER_OF_MAXIMUM_TRIGGERS> tauTriggerBitSet_woL3_2_;

      HLTConfigProvider hltConfig_;

      std::vector <std::string> triggerModules_;
      TString filterLabel_;
      TTree* filterLabelsTree_;
      unsigned int lastFilterInd_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
NtuplizerVBF::NtuplizerVBF(const edm::ParameterSet& iConfig) :
    tauTag_ (consumes<pat::TauRefVector>  (iConfig.getParameter<edm::InputTag>("taus"))),
    leadJet_ (consumes<pat::JetRefVector>  (iConfig.getParameter<edm::InputTag>("jets"))),
    jetIdTag_ (consumes<edm::ValueMap<int>> (iConfig.getParameter<edm::InputTag>("jetId"))),
    genParticlesTag_ (consumes<edm::View<pat::GenericParticle>>  (iConfig.getParameter<edm::InputTag>("genParticles"))),
    triggerObjects_ (consumes<pat::TriggerObjectStandAloneCollection> (iConfig.getParameter<edm::InputTag>("triggerSet"))),
    triggerBits_ (consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("triggerResultsLabel"))),
    VtxTag_ (consumes<std::vector<reco::Vertex>> (iConfig.getParameter<edm::InputTag>("Vertices"))),
    puTag_ (consumes<std::vector<PileupSummaryInfo>> (iConfig.getParameter<edm::InputTag>("puInfo"))),
    processName_ (iConfig.getParameter<edm::InputTag>("triggerResultsLabel")),
    treeName_ (iConfig.getParameter<std::string>("treeName")),
    isMC_ (iConfig.getParameter<bool>("isMC")),
    filterPath_ (iConfig.getParameter<std::string>("filterPath"))
{
   TString triggerName;
   edm::Service<TFileService> fs;
   triggerNamesTree_ = fs->make<TTree>("triggerNames", "triggerNames");
   triggerNamesTree_->Branch("triggerNames", &triggerName);

   filterLabelsTree_ = fs -> make<TTree>("filterLabels", "filterLabels");
   filterLabelsTree_ -> Branch("filterLabels", &filterLabel_);


   //Building the trigger arrays
   const std::vector<edm::ParameterSet>& HLTList = iConfig.getParameter <std::vector<edm::ParameterSet> > ("triggerListProbe");
   for (const edm::ParameterSet& parameterSet : HLTList) {
       tParameterSet pSet;
       pSet.hltPath = parameterSet.getParameter<std::string>("HLT");
       triggerName = pSet.hltPath;
       pSet.hltFilters1 = parameterSet.getParameter<std::vector<std::string> >("path1");
       pSet.hltFilters2 = parameterSet.getParameter<std::vector<std::string> >("path2");
       pSet.leg1 = parameterSet.getParameter<int>("leg1");
       pSet.leg2 = parameterSet.getParameter<int>("leg2");
       parameters_.push_back(pSet);

       triggerNamesTree_->Fill();
   }

   Initialize();
   return;

}


NtuplizerVBF::~NtuplizerVBF()
{
    return;
}


//
// member functions
//

// ------------ method called for each event  ------------
void NtuplizerVBF::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    Initialize();

    indexevents_ = iEvent.id().event();
    runNumber_ = iEvent.id().run();
    lumi_ = iEvent.luminosityBlock();

    edm::Handle<pat::TauRefVector> taus;
    edm::Handle<pat::JetRefVector> jets;
    edm::Handle<edm::ValueMap<int>> jet_id_decisions; 
    edm::Handle<edm::View<pat::GenericParticle> > genParticles;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<std::vector<reco::Vertex> >  vertices;
    edm::Handle<std::vector<PileupSummaryInfo>> puInfo;

    iEvent.getByToken(tauTag_, taus);
    iEvent.getByToken(leadJet_, jets);
    iEvent.getByToken(jetIdTag_, jet_id_decisions);
    iEvent.getByToken(triggerObjects_, triggerObjects);
    iEvent.getByToken(triggerBits_, triggerBits);
    iEvent.getByToken(VtxTag_, vertices);
    iEvent.getByToken(puTag_, puInfo);

    if(isMC_)
      iEvent.getByToken(genParticlesTag_, genParticles);

    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

    const pat::JetRef leadJet = jets->at(0);
    leadJetPt_ = leadJet->pt();
    leadJetEta_ = leadJet->eta();
    leadJetPhi_ = leadJet->phi();

    const pat::JetRef trailJet = jets->at(1);
    trailJetPt_ = trailJet->pt();
    trailJetEta_ = trailJet->eta();
    trailJetPhi_ = trailJet->phi();

    if (jets->size() > 2)
    {
        const pat::JetRef addJet = jets->at(2);
        thirdJetPt_ = addJet->pt();
        thirdJetEta_ = addJet->eta();
        thirdJetPhi_ = addJet->phi();
    }

    Mjj_ = (leadJet->p4() + trailJet->p4()).M();
    nJets_ = jets->size();

    Nvtx_ = vertices->size();

    nTruePU_ = -99;
    if (isMC_)
    {
        std::vector<PileupSummaryInfo>::const_iterator PVI;
        for(PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI)
        {
            if(PVI->getBunchCrossing() == 0)
            {
                float nTrueInt = PVI->getTrueNumInteractions();
                nTruePU_ = nTrueInt;
                break;
            }
        }
    }

    assert(taus->size() == 2);
    const pat::TauRef tauLead = taus->at(0);
    const pat::TauRef tauTrail = taus->at(1);

    tauTriggerBitSet_1_.reset();
    tauTriggerBitSet_woL3_1_.reset();
    tauTriggerBitSet_2_.reset();
    tauTriggerBitSet_woL3_2_.reset();

    bool leadJetL1Matched = false;
    bool trailJetL1Matched = false;

    for (const pat::TriggerObjectStandAlone obj: *triggerObjects)
    {
        // TODO update filter matching
        const float dR1 = deltaR(*tauLead, obj);
        if (dR1 < 0.5)
        {
            const edm::TriggerNames::Strings& triggerNames = names.triggerNames();
            // Looking for the path index.
            unsigned int x = 0;
            unsigned int y = 0;
            for (const tParameterSet& parameter : parameters_)
            {
                if ((parameter.hltPathIndex >= 0) && obj.hasPathName(triggerNames[parameter.hltPathIndex], false, true))
                {
                    const std::vector<std::string>& filters = parameter.hltFilters1;
                    if (hasFilters(obj, filters))
                    {
                        tauTriggerBitSet_1_[x] = true;
                    }
                }
                x += 1;
                if ((parameter.hltPathIndex >= 0) && obj.hasPathName(triggerNames[parameter.hltPathIndex], false, false))
                {
                    const std::vector<std::string>& filters = parameter.hltFilters1;
                    if (hasFilters(obj, filters))
                    {
                        tauTriggerBitSet_woL3_1_[y] = true;
                    }
                }
                y += 1;
            }
        }
        const float dR2 = deltaR(*tauTrail, obj);
        if (dR2 < 0.5)
        {
            const edm::TriggerNames::Strings& triggerNames = names.triggerNames();
            // Looking for the path index.
            unsigned int x = 0;
            unsigned int y = 0;
            for (const tParameterSet& parameter : parameters_)
            {
                if ((parameter.hltPathIndex >= 0) && obj.hasPathName(triggerNames[parameter.hltPathIndex], false, true))
                {
                    const std::vector<std::string>& filters = parameter.hltFilters1;
                    if (hasFilters(obj, filters))
                    {
                        tauTriggerBitSet_2_[x] = true;
                    }
                }
                x += 1;
                if ((parameter.hltPathIndex >= 0) && obj.hasPathName(triggerNames[parameter.hltPathIndex], false, false))
                {
                    const std::vector<std::string>& filters = parameter.hltFilters1;
                    if (hasFilters(obj, filters))
                    {
                        tauTriggerBitSet_woL3_2_[y] = true;
                    }
                }
                y += 1;
            }
        }
        std::vector<std::string> filters = {"hltL1VBFDiJetOR"};
        const float dRJet1 = deltaR(*leadJet, obj);
        if (dRJet1 < 0.5)
        {
            if (obj.hasPathName(filterPath_+"*", false, false))
            {
                if (hasFilters(obj, filters))
                {
                    leadJetL1Matched = true;
                }
            }
        }
        const float dRJet2 = deltaR(*trailJet, obj);
        if (dRJet2 < 0.5)
        {
            if (obj.hasPathName(filterPath_+"*", false, false))
            {
                if (hasFilters(obj, filters))
                {
                    trailJetL1Matched = true;
                }
            }
        }
        // std::cout << "Kinematics: " << "pT: " << obj.pt() << " eta: " << obj.eta() << " phi: " << obj.phi() << std::endl;
        // bool isBoth = obj.hasPathName( filterPath_ + "*", true, true );
        // bool isL3   = obj.hasPathName( filterPath_ + "*", false, true );
        // bool isLF   = obj.hasPathName( filterPath_ + "*", true, false );
        // bool isNone = obj.hasPathName( filterPath_ + "*", false, false );
        // std::cout << "   " << filterPath_ + "*";
        // if (isBoth) std::cout << "(L,3)";
        // if (isL3 && !isBoth) std::cout << "(*,3)";
        // if (isLF && !isBoth) std::cout << "(L,*)";
        // if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
        //std::cout << std::endl;
        if (obj.hasPathName(filterPath_ + "*", false, false))
        {
            for (std::vector<std::string>::reverse_iterator filterName = triggerModules_.rbegin(); filterName != triggerModules_.rend(); filterName+=1)
            {
                if (obj.hasFilterLabel(*filterName))
                {
                    if (triggerModules_.rend() - filterName > lastFilter_)
                    {                                                                                                                                                                                               
                            lastFilter_ = triggerModules_.rend() - filterName;
                    }
                }
            }
        }
    }   

    tauTriggerBits_1_ = tauTriggerBitSet_1_.to_ulong();
    tauTriggerBits_woL3_1_ = tauTriggerBitSet_woL3_1_.to_ulong();
    tauTriggerBits_2_ = tauTriggerBitSet_2_.to_ulong();
    tauTriggerBits_woL3_2_ = tauTriggerBitSet_woL3_2_.to_ulong();

    jetPairL1Matched_ = leadJetL1Matched && trailJetL1Matched;

    tauPt_1_ = tauLead->pt();
    tauEta_1_ = tauLead->eta();
    tauPhi_1_ = tauLead->phi();
    tauCharge_1_ = tauLead->charge();
    tauDM_1_ = tauLead->decayMode();
    tauTrkPt_1_ = tauLead->leadChargedHadrCand()->pt();
    tauPt_2_ = tauTrail->pt();
    tauEta_2_ = tauTrail->eta();
    tauPhi_2_ = tauTrail->phi();
    tauCharge_2_ = tauTrail->charge();
    tauDM_2_ = tauTrail->decayMode();
    tauTrkPt_2_ = tauTrail->leadChargedHadrCand()->pt();
    if (isMC_)
    {
        const edm::View<pat::GenericParticle>* genparts = genParticles.product();
        tau_genindex_1_ = GenIndex(tauLead, genparts);
        tau_genindex_2_ = GenIndex(tauTrail, genparts);
    }
    tauByVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1_ = tauLead->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017");
    tauByVLooseIsolationMVArun2017v2DBoldDMwLT2017_1_ = tauLead->tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017");
    tauByLooseIsolationMVArun2017v2DBoldDMwLT2017_1_ = tauLead->tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017");
    tauByMediumIsolationMVArun2017v2DBoldDMwLT2017_1_ = tauLead->tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017");
    tauByTightIsolationMVArun2017v2DBoldDMwLT2017_1_ = tauLead->tauID("byTightIsolationMVArun2017v2DBoldDMwLT2017");
    tauByVTightIsolationMVArun2017v2DBoldDMwLT2017_1_ = tauLead->tauID("byVTightIsolationMVArun2017v2DBoldDMwLT2017");
    tauByVVTightIsolationMVArun2017v2DBoldDMwLT2017_1_ = tauLead->tauID("byVVTightIsolationMVArun2017v2DBoldDMwLT2017");
    tauAgainstMuonLoose3_1_ = tauLead->tauID("againstMuonLoose3");
    tauAgainstMuonTight3_1_ = tauLead->tauID("againstMuonTight3");
    tauAgainstElectronVLooseMVA6_1_ = tauLead->tauID("againstElectronVLooseMVA6");
    tauAgainstElectronLooseMVA6_1_ = tauLead->tauID("againstElectronLooseMVA6");
    tauAgainstElectronMediumMVA6_1_ = tauLead->tauID("againstElectronMediumMVA6");
    tauAgainstElectronTightMVA6_1_ = tauLead->tauID("againstElectronTightMVA6");
    tauAgainstElectronVTightMVA6_1_ = tauLead->tauID("againstElectronVTightMVA6");

    tauByVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2_ = tauTrail->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017");
    tauByVLooseIsolationMVArun2017v2DBoldDMwLT2017_2_ = tauTrail->tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017");
    tauByLooseIsolationMVArun2017v2DBoldDMwLT2017_2_ = tauTrail->tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017");
    tauByMediumIsolationMVArun2017v2DBoldDMwLT2017_2_ = tauTrail->tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017");
    tauByTightIsolationMVArun2017v2DBoldDMwLT2017_2_ = tauTrail->tauID("byTightIsolationMVArun2017v2DBoldDMwLT2017");
    tauByVTightIsolationMVArun2017v2DBoldDMwLT2017_2_ = tauTrail->tauID("byVTightIsolationMVArun2017v2DBoldDMwLT2017");
    tauByVVTightIsolationMVArun2017v2DBoldDMwLT2017_2_ = tauTrail->tauID("byVVTightIsolationMVArun2017v2DBoldDMwLT2017");
    tauAgainstMuonLoose3_2_ = tauTrail->tauID("againstMuonLoose3");
    tauAgainstMuonTight3_2_ = tauTrail->tauID("againstMuonTight3");
    tauAgainstElectronVLooseMVA6_2_ = tauTrail->tauID("againstElectronVLooseMVA6");
    tauAgainstElectronLooseMVA6_2_ = tauTrail->tauID("againstElectronLooseMVA6");
    tauAgainstElectronMediumMVA6_2_ = tauTrail->tauID("againstElectronMediumMVA6");
    tauAgainstElectronTightMVA6_2_ = tauTrail->tauID("againstElectronTightMVA6");
    tauAgainstElectronVTightMVA6_2_ = tauTrail->tauID("againstElectronVTightMVA6");

    //std::cout << "Fill event with: EventNumber " << indexevents_ << " RunNumber " << runNumber_ << " and LumiSection: " << lumi_ << std::endl;
    tree_->Fill();
    // else
    // {
    //     std::cout << "Did not fill event with: EventNumber " << indexevents_ << " RunNumber " << runNumber_ << " and LumiSection: " << lumi_ << std::endl;
    //     std::cout << "reason: " << isTagHLTmatched ? std::cout << " no tau found" : std::cout << " tag not HLT matched";
    //     std::cout << std::endl;
    // }
}


// ------------ method called once each job just before starting event loop  ------------
void NtuplizerVBF::beginJob()
{
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>(treeName_.c_str(), treeName_.c_str());

    tree_->Branch("EventNumber", &indexevents_,"EventNumber/l");
    tree_->Branch("RunNumber", &runNumber_, "RunNumber/I");
    tree_->Branch("lumi", &lumi_, "lumi/I");

    tree_->Branch("tauPt_1", &tauPt_1_, "tauPt_1/F");
    tree_->Branch("tauEta_1", &tauEta_1_, "tauEta_1/F");
    tree_->Branch("tauPhi_1", &tauPhi_1_, "tauPhi_1/F");
    tree_->Branch("tauCharge_1", &tauCharge_1_, "tauCharge_1/I");
    tree_->Branch("tauDM_1", &tauDM_1_, "tauDM_1/I");
    tree_->Branch("tauTrkPt_1", &tauTrkPt_1_, "tauTrkPt_1/F");
    tree_->Branch("tau_genindex_1", &tau_genindex_1_, "tau_genindex_1/I");

    tree_->Branch("tauPt_2", &tauPt_2_, "tauPt_2/F");
    tree_->Branch("tauEta_2", &tauEta_2_, "tauEta_2/F");
    tree_->Branch("tauPhi_2", &tauPhi_2_, "tauPhi_2/F");
    tree_->Branch("tauCharge_2", &tauCharge_2_, "tauCharge_2/I");
    tree_->Branch("tauDM_2", &tauDM_2_, "tauDM_2/I");
    tree_->Branch("tauTrkPt_2", &tauTrkPt_2_, "tauTrkPt_2/F");
    tree_->Branch("tau_genindex_2", &tau_genindex_2_, "tau_genindex_2/I");

    tree_->Branch("tauByVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1", &tauByVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1_, "tauByVVLooseIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("tauByVLooseIsolationMVArun2017v2DBoldDMwLT2017_1", &tauByVLooseIsolationMVArun2017v2DBoldDMwLT2017_1_, "tauByVLooseIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("tauByLooseIsolationMVArun2017v2DBoldDMwLT2017_1", &tauByLooseIsolationMVArun2017v2DBoldDMwLT2017_1_, "tauByLooseIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("tauByMediumIsolationMVArun2017v2DBoldDMwLT2017_1", &tauByMediumIsolationMVArun2017v2DBoldDMwLT2017_1_, "tauByMediumIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("tauByTightIsolationMVArun2017v2DBoldDMwLT2017_1", &tauByTightIsolationMVArun2017v2DBoldDMwLT2017_1_, "tauByTightIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("tauByVTightIsolationMVArun2017v2DBoldDMwLT2017_1", &tauByVTightIsolationMVArun2017v2DBoldDMwLT2017_1_, "tauByVTightIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("tauByVVTightIsolationMVArun2017v2DBoldDMwLT2017_1", &tauByVVTightIsolationMVArun2017v2DBoldDMwLT2017_1_, "tauByVVTightIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("tauAgainstMuonLoose3_1", &tauAgainstMuonLoose3_1_,"tauAgainstMuonLoose3_1/O");
    tree_->Branch("tauAgainstMuonTight3_1", &tauAgainstMuonTight3_1_, "tauAgainstMuonTight3_1/O");
    tree_->Branch("tauAgainstElectronVLooseMVA6_1", &tauAgainstElectronVLooseMVA6_1_, "tauAgainstElectronVLooseMVA6_1/O");
    tree_->Branch("tauAgainstElectronLooseMVA6_1", &tauAgainstElectronLooseMVA6_1_, "tauAgainstElectronLooseMVA6_1/O");
    tree_->Branch("tauAgainstElectronMediumMVA6_1", &tauAgainstElectronMediumMVA6_1_, "tauAgainstElectronMediumMVA6_1/O");
    tree_->Branch("tauAgainstElectronTightMVA6_1", &tauAgainstElectronTightMVA6_1_, "tauAgainstElectronTightMVA6_1/O");
    tree_->Branch("tauAgainstElectronVTightMVA6_1", &tauAgainstElectronVTightMVA6_1_, "tauAgainstElectronVTightMVA6_1/O");
    tree_->Branch("tauByVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2", &tauByVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2_, "tauByVVLooseIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("tauByVLooseIsolationMVArun2017v2DBoldDMwLT2017_2", &tauByVLooseIsolationMVArun2017v2DBoldDMwLT2017_2_, "tauByVLooseIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("tauByLooseIsolationMVArun2017v2DBoldDMwLT2017_2", &tauByLooseIsolationMVArun2017v2DBoldDMwLT2017_2_, "tauByLooseIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("tauByMediumIsolationMVArun2017v2DBoldDMwLT2017_2", &tauByMediumIsolationMVArun2017v2DBoldDMwLT2017_2_, "tauByMediumIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("tauByTightIsolationMVArun2017v2DBoldDMwLT2017_2", &tauByTightIsolationMVArun2017v2DBoldDMwLT2017_2_, "tauByTightIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("tauByVTightIsolationMVArun2017v2DBoldDMwLT2017_2", &tauByVTightIsolationMVArun2017v2DBoldDMwLT2017_2_, "tauByVTightIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("tauByVVTightIsolationMVArun2017v2DBoldDMwLT2017_2", &tauByVVTightIsolationMVArun2017v2DBoldDMwLT2017_2_, "tauByVVTightIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("tauAgainstMuonLoose3_2", &tauAgainstMuonLoose3_2_,"tauAgainstMuonLoose3_2/O");
    tree_->Branch("tauAgainstMuonTight3_2", &tauAgainstMuonTight3_2_, "tauAgainstMuonTight3_2/O");
    tree_->Branch("tauAgainstElectronVLooseMVA6_2", &tauAgainstElectronVLooseMVA6_2_, "tauAgainstElectronVLooseMVA6_2/O");
    tree_->Branch("tauAgainstElectronLooseMVA6_2", &tauAgainstElectronLooseMVA6_2_, "tauAgainstElectronLooseMVA6_2/O");
    tree_->Branch("tauAgainstElectronMediumMVA6_2", &tauAgainstElectronMediumMVA6_2_, "tauAgainstElectronMediumMVA6_2/O");
    tree_->Branch("tauAgainstElectronTightMVA6_2", &tauAgainstElectronTightMVA6_2_, "tauAgainstElectronTightMVA6_2/O");
    tree_->Branch("tauAgainstElectronVTightMVA6_2", &tauAgainstElectronVTightMVA6_2_, "tauAgainstElectronVTightMVA6_2/O");
    tree_->Branch("tauTriggerBits_1", &tauTriggerBits_1_, "tauTriggerBits_1/l");
    tree_->Branch("tauTriggerBits_woL3_1", &tauTriggerBits_woL3_1_, "tauTriggerBits_woL3_1/l");
    tree_->Branch("tauTriggerBits_2", &tauTriggerBits_2_, "tauTriggerBits_2/l");
    tree_->Branch("tauTriggerBits_woL3_2", &tauTriggerBits_woL3_2_, "tauTriggerBits_woL3_2/l");

    tree_->Branch("leadJetPt", &leadJetPt_, "leadJetPt/F");
    tree_->Branch("leadJetEta", &leadJetEta_, "leadJetEta/F");
    tree_->Branch("leadJetPhi", &leadJetPhi_, "leadJetPhi/F");

    tree_->Branch("trailJetPt", &trailJetPt_, "trailJetPt/F");
    tree_->Branch("trailJetEta", &trailJetEta_, "trailJetEta/F");
    tree_->Branch("trailJetPhi", &trailJetPhi_, "trailJetPhi/F");

    tree_->Branch("thirdJetPt", &thirdJetPt_, "thirdJetPt/F");
    tree_->Branch("thirdJetEta", &thirdJetEta_, "thirdJetEta/F");
    tree_->Branch("thirdJetPhi", &thirdJetPhi_, "thirdJetPhi/F");

    tree_->Branch("jetPairL1Matched", &jetPairL1Matched_, "jetPairL1Matched/O");
    
    tree_->Branch("nJets", &nJets_, "nJets/I");
    tree_->Branch("Mjj", &Mjj_, "Mjj/F");

    tree_->Branch("Nvtx", &Nvtx_, "Nvtx/I");
    tree_->Branch("nTruePU", &nTruePU_, "nTruePU/F");

    tree_ -> Branch("lastFilter", &lastFilter_, "lastFilter/I");

    return;
}

// ------------ method called once each job just after ending the event loop  ------------
void NtuplizerVBF::endJob() 
{
    return;
}


void NtuplizerVBF::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    Bool_t changedConfig = false;

    if (!hltConfig_.init(iRun, iSetup, processName_.process(), changedConfig))
    {
        edm::LogError("HLTMatchingFilter") << "Initialization of HLTConfigProvider failed!!";
        return;
    }

    const edm::TriggerNames::Strings& triggerNames = hltConfig_.triggerNames();
    std::cout << " ===== LOOKING FOR THE PATH INDEXES =====" << std::endl;

    for (tParameterSet& parameter : parameters_){
        const std::string& hltPath = parameter.hltPath;
        bool found = false;
        for(unsigned int j=0; j < triggerNames.size(); j++)
        {
            if (triggerNames[j].find(hltPath) != std::string::npos) {
                found = true;
                parameter.hltPathIndex = j;

                std::cout << "### FOUND AT INDEX #" << j << " --> " << triggerNames[j] << std::endl;
                // Look for the trigger filters running in this configuration.
                if (hltPath==filterPath_)
                {
                    lastFilterInd_ = j;
                }
            }
        }
        if (!found) parameter.hltPathIndex = -1;
    }
    // Get trigger modules which ran with saveTags option, e.g. important EDFilters                                                                                                                                 
     triggerModules_ = hltConfig_.saveTagsModules(lastFilterInd_);
     for (const std::string triggerModule: triggerModules_)
     {
         filterLabel_ = triggerModule;
         filterLabelsTree_->Fill();
     }

    return;
}

void NtuplizerVBF::Initialize()
{
    indexevents_ = 0;
    runNumber_ = 0;
    lumi_ = 0;

    tauPt_1_ = -1.;
    tauEta_1_ = -999;
    tauPhi_1_ = -999;
    tauCharge_1_ = 0;
    tauDM_1_ = -1;
    tauTrkPt_1_ = -1.;
    tau_genindex_1_ = -1;
    tauPt_2_ = -1.;
    tauEta_2_ = -999;
    tauPhi_2_ = -999;
    tauCharge_2_ = 0;
    tauDM_2_ = -1;
    tauTrkPt_2_ = -1.;
    tau_genindex_2_ = -1;

    tauByVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1_ = 0;
    tauByVLooseIsolationMVArun2017v2DBoldDMwLT2017_1_ = 0;
    tauByLooseIsolationMVArun2017v2DBoldDMwLT2017_1_ = 0;
    tauByMediumIsolationMVArun2017v2DBoldDMwLT2017_1_ = 0;
    tauByTightIsolationMVArun2017v2DBoldDMwLT2017_1_ = 0;
    tauByVTightIsolationMVArun2017v2DBoldDMwLT2017_1_ = 0;
    tauByVVTightIsolationMVArun2017v2DBoldDMwLT2017_1_ = 0;
    tauAgainstMuonLoose3_1_ = 0;
    tauAgainstMuonTight3_1_ = 0;
    tauAgainstElectronVLooseMVA6_1_ = 0;
    tauAgainstElectronLooseMVA6_1_ = 0;
    tauAgainstElectronMediumMVA6_1_ = 0;
    tauAgainstElectronTightMVA6_1_ = 0;
    tauAgainstElectronVTightMVA6_1_ = 0;
    tauByVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2_ = 0;
    tauByVLooseIsolationMVArun2017v2DBoldDMwLT2017_2_ = 0;
    tauByLooseIsolationMVArun2017v2DBoldDMwLT2017_2_ = 0;
    tauByMediumIsolationMVArun2017v2DBoldDMwLT2017_2_ = 0;
    tauByTightIsolationMVArun2017v2DBoldDMwLT2017_2_ = 0;
    tauByVTightIsolationMVArun2017v2DBoldDMwLT2017_2_ = 0;
    tauByVVTightIsolationMVArun2017v2DBoldDMwLT2017_2_ = 0;
    tauAgainstMuonLoose3_2_ = 0;
    tauAgainstMuonTight3_2_ = 0;
    tauAgainstElectronVLooseMVA6_2_ = 0;
    tauAgainstElectronLooseMVA6_2_ = 0;
    tauAgainstElectronMediumMVA6_2_ = 0;
    tauAgainstElectronTightMVA6_2_ = 0;
    tauAgainstElectronVTightMVA6_2_ = 0;


    jetPairL1Matched_ = false;
    leadJetPt_ = -1.;
    leadJetEta_ = -999;
    leadJetPhi_ = -999;

    trailJetPt_ = -1.;
    trailJetEta_ = -999;
    trailJetPhi_ = -999;

    thirdJetPt_ = -1.;
    thirdJetEta_ = -999;
    thirdJetPhi_ = -999;
    Mjj_ = -1.;
    nJets_ = -1;


    Nvtx_ = 0;
    nTruePU_ = 0.;

    tauTriggerBits_1_ = 999;
    tauTriggerBits_woL3_1_ = 999;
    tauTriggerBits_2_ = 999;
    tauTriggerBits_woL3_2_ = 999;

    lastFilter_ = 0;
    return; 
}

void NtuplizerVBF::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    return;
}

bool NtuplizerVBF::hasFilters(const pat::TriggerObjectStandAlone&  obj , const std::vector<std::string>& filtersToLookFor)
{
    const std::vector<std::string>& eventLabels = obj.filterLabels();
    for (const std::string& filter : filtersToLookFor)
    {
        //Looking for matching filters
        bool found = false;
        for (const std::string& label : eventLabels)
        {
            //if (label == std::string("hltOverlapFilterIsoMu17MediumIsoPFTau40Reg"))
            if (label == filter)
            {
                //std::cout << "#### FOUND FILTER " << label << " == " << filter << " ####" << std::endl;
                found = true;
            }
        }
        if(!found)
        {
            return false;
        }
    }
    return true;
}

int NtuplizerVBF::GenIndex(const pat::TauRef& tau, const edm::View<pat::GenericParticle>* genparts)
{
    float dRmin = 1.0;
    int genMatchInd = -1;

    for(edm::View<pat::GenericParticle>::const_iterator genpart = genparts->begin(); genpart!=genparts->end();++genpart)
    {

        int flags = genpart->userInt ("generalGenFlags");
        int apdg = abs(genpart->pdgId());
        float pT = genpart->p4().pt();

        if( !( apdg==11 || apdg==13 || apdg==66615) ) continue;

        if( apdg==11 || apdg==13)
        {
            if( !(pT>8 && (flags&1 || (flags>>5)&1)) ) continue;
        }
        else if(apdg==66615)
        {
            int tauMothInd = genpart->userInt("TauMothIndex");
            pat::GenericParticle mother = (*genparts)[tauMothInd];
            int flags_tau = mother.userInt ("generalGenFlags");
            if( !(pT>15 && flags_tau&1) ) continue;
        }

        float dR = deltaR(*tau,*genpart);
        if(dR<0.2 && dR<dRmin)
        {
            dRmin = dR;
            if(apdg==11)
            {
                if(flags&1) genMatchInd = 1;
                else if((flags>>5)&1) genMatchInd = 3;
            }
            else if(apdg==13)
            {
                if(flags&1) genMatchInd = 2;
                else if((flags>>5)&1) genMatchInd = 4;
            }
            else if(apdg==66615)
                genMatchInd = 5;
        }

    }

    return genMatchInd;
}

//define this as a plug-in
#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(NtuplizerVBF);
#endif
