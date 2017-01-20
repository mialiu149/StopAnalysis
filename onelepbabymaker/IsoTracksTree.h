#ifndef ISOTRACKSTREE_H
#define ISOTRACKSTREE_H

#include <vector>
#include <string>
#include "Math/LorentzVector.h"

class TTree;

// typedefs
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef std::vector<LorentzVector> vecLorentzVector;
typedef std::vector<float> vecd;
typedef std::vector<std::string> vecs;
typedef std::vector<int> veci;
typedef std::vector<bool> vecb;

class IsoTracksTree
{
public:
    IsoTracksTree ();
    virtual ~IsoTracksTree () {}

    void Reset ();
    void SetBranches (TTree* tree);
    void SetAliases (TTree* tree) const;
    void FillCommon (int idx);
    int  nhighPtPFcands;
    vecLorentzVector isoTracks_p4;
    veci isoTracks_charge;
    vecd isoTracks_absIso;
    vecd isoTracks_dz;
    veci isoTracks_pdgId;
    vecb isoTracks_isVetoTrack;
    vecb isoTracks_isVetoTrack_v2;
    vecb isoTracks_isVetoTrack_v3;
    vecd highPtPFcands_pt;
    vecd highPtPFcands_eta;
    vecd highPtPFcands_phi;
    vecd highPtPFcands_mass;
    vecd highPtPFcands_absIso;
    vecd highPtPFcands_relIsoAn04;
    vecd highPtPFcands_dz;
    vecd highPtPFcands_pdgId;
    vecd highPtPFcands_mcMatchId;
};

#endif
