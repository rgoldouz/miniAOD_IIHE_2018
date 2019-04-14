#ifndef ETSORT
#define ETSORT

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
inline bool gsfEtGreater(const reco::GsfElectron &gsf1,const reco::GsfElectron &gsf2){
  float et1 = gsf1.caloEnergy() * sin(gsf1.p4().theta());
  float et2 = gsf2.caloEnergy() * sin(gsf2.p4().theta());
  return (et1 > et2);
}

inline bool scEtGreater(const reco::SuperCluster *sc1,const reco::SuperCluster *sc2) {
  return (((sc1->energy() + sc1->preshowerEnergy()) )/cosh((sc1)->eta()) >((sc2->energy() + sc2->preshowerEnergy()) )/cosh((sc2)->eta()) );
}

inline bool refScEtGreater(reco::SuperClusterRef sc1,reco::SuperClusterRef sc2) {
  return (((sc1->energy() + sc1->preshowerEnergy()) )/cosh((sc1)->eta()) >((sc2->energy() + sc2->preshowerEnergy()) )/cosh((sc2)->eta()));
}

#endif

