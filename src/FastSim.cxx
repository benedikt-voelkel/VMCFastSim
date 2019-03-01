// Some license stuff

#include <TError.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <TVirtualMCApplication.h>
#include <TVirtualMCStack.h>
#include <TMCManagerStack.h>
#include <TParticle.h>
#include <TMCParticleStatus.h>

#include "VMCFastSim/FastSim.h"

using namespace vmcfastsim::base;


template <typename K>
Bool_t FastSim<K>::IsRootGeometrySupported() const
{
  return fIsRootGeometrySupported;
}

template <typename K>
void FastSim<K>::Material(Int_t& kmat, const char* name, Double_t a, Double_t z,
                       Double_t dens, Double_t radl, Double_t absl,
                       Float_t* buf, Int_t nwbuf)
{
  Warning("Material", "Not supported");
}

template <typename K>
void FastSim<K>::Material(Int_t& kmat, const char* name, Double_t a,
                       Double_t z, Double_t dens, Double_t radl, Double_t absl,
                       Double_t* buf, Int_t nwbuf)
{
  Warning("Material", "Not supported");
}

template <typename K>
void FastSim<K>::Mixture(Int_t& kmat, const char *name, Float_t *a, Float_t *z,
                      Double_t dens, Int_t nlmat, Float_t *wmat)
{
  Warning("Mixture", "Not supported");
}

template <typename K>
void FastSim<K>::Mixture(Int_t& kmat, const char *name, Double_t *a,
                      Double_t *z, Double_t dens, Int_t nlmat, Double_t *wmat)
{
  Warning("Mixture", "Not supported");
}

template <typename K>
void FastSim<K>::Medium(Int_t& kmed, const char *name, Int_t nmat,
                     Int_t isvol, Int_t ifield, Double_t fieldm,
                     Double_t tmaxfd, Double_t stemax, Double_t deemax,
                     Double_t epsil, Double_t stmin, Float_t* ubuf, Int_t nbuf)
{
  Warning("Medium", "Not supported");
}


template <typename K>
void FastSim<K>::Medium(Int_t& kmed, const char *name, Int_t nmat, Int_t isvol,
                     Int_t ifield, Double_t fieldm, Double_t tmaxfd,
                     Double_t stemax, Double_t deemax, Double_t epsil,
                     Double_t stmin, Double_t* ubuf, Int_t nbuf)
{
  Warning("Medium", "Not supported");
}

template <typename K>
void FastSim<K>::Matrix(Int_t& krot, Double_t thetaX, Double_t phiX,
                     Double_t thetaY, Double_t phiY, Double_t thetaZ,
                     Double_t phiZ)
{
  Warning("Matrix", "Not supported");
}

template <typename K>
void FastSim<K>::Gstpar(Int_t itmed, const char *param, Double_t parval)
{
  Warning("Gstpar", "Not yet implemented");
}

template <typename K>
Int_t FastSim<K>::Gsvolu(const char *name, const char *shape, Int_t nmed,
                      Float_t *upar, Int_t np)
{
  Warning("Gsvolu", "Not supported");
  return -1;
}

template <typename K>
Int_t FastSim<K>::Gsvolu(const char *name, const char *shape, Int_t nmed,
                      Double_t *upar, Int_t np)
{
  Warning("Gsvolu", "Not supported");
  return -1;
}

template <typename K>
void FastSim<K>::Gsdvn(const char *name, const char *mother, Int_t ndiv,
                    Int_t iaxis)
{
  Warning("Gsdvn", "Not supported");
}


template <typename K>
void FastSim<K>::Gsdvn2(const char *name, const char *mother, Int_t ndiv,
                     Int_t iaxis, Double_t c0i, Int_t numed)
{
  Warning("Gsdvn2", "Not supported");
}


template <typename K>
void FastSim<K>::Gsdvt(const char *name, const char *mother, Double_t step,
                    Int_t iaxis, Int_t numed, Int_t ndvmx)
{
  Warning("Gsdvt", "Not supported");
}


template <typename K>
void FastSim<K>::Gsdvt2(const char *name, const char *mother, Double_t step,
                     Int_t iaxis, Double_t c0, Int_t numed, Int_t ndvmx)
{
  Warning("Gsdvt2", "Not supported");
}

template <typename K>
void FastSim<K>::Gsord(const char *name, Int_t iax)
{
  Warning("Gsord", "Not supported");
}

template <typename K>
void FastSim<K>::Gspos(const char *name, Int_t nr, const char *mother,
                    Double_t x, Double_t y, Double_t z, Int_t irot,
                    const char *konly)
{
  Warning("Gspos", "Not supported");
}

template <typename K>
void FastSim<K>::Gsposp(const char *name, Int_t nr, const char *mother,
                     Double_t x, Double_t y, Double_t z, Int_t irot,
                     const char *konly, Float_t *upar, Int_t np)
{
  Warning("Gsposp", "Not supported");
}

template <typename K>
void FastSim<K>::Gsposp(const char *name, Int_t nr, const char *mother,
                     Double_t x, Double_t y, Double_t z, Int_t irot,
                     const char *konly, Double_t *upar, Int_t np)
{
  Warning("Gsposp", "Not supported");
}

template <typename K>
void FastSim<K>::Gsbool(const char* onlyVolName, const char* manyVolName)
{
  Warning("Gsbool", "Not supported");
}

template <typename K>
void FastSim<K>::SetCerenkov(Int_t itmed, Int_t npckov, Float_t *ppckov,
                          Float_t *absco, Float_t *effic, Float_t *rindex)
{
  Warning("Gsbool", "Not supported");
}

template <typename K>
void FastSim<K>::SetCerenkov(Int_t itmed, Int_t npckov, Double_t *ppckov,
                          Double_t *absco, Double_t *effic, Double_t *rindex)
{
  Warning("Gsbool", "Not supported");
}

template <typename K>
void FastSim<K>::DefineOpSurface(const char* name,
                              EMCOpSurfaceModel model,
                              EMCOpSurfaceType surfaceType,
                              EMCOpSurfaceFinish surfaceFinish,
                              Double_t sigmaAlpha)
{
  Warning("DefineOpSurface", "Not supported");
}

template <typename K>
void FastSim<K>::SetBorderSurface(const char* name,
                               const char* vol1Name, int vol1CopyNo,
                               const char* vol2Name, int vol2CopyNo,
                               const char* opSurfaceName)
{
  Warning("SetBorderSurface", "Not supported");
}

template <typename K>
void FastSim<K>::SetSkinSurface(const char* name,
                             const char* volName,
                             const char* opSurfaceName)
{
  Warning("SetSkinSurface", "Not supported");
}

template <typename K>
void FastSim<K>::SetMaterialProperty(Int_t itmed, const char* propertyName,
                                  Int_t np, Double_t* pp, Double_t* values)
{
  Warning("SetMaterialProperty", "Not supported");
}

template <typename K>
void FastSim<K>::SetMaterialProperty(Int_t itmed, const char* propertyName,
                                  Double_t value)
{
  Warning("SetMaterialProperty", "Not supported");
}

template <typename K>
void FastSim<K>::SetMaterialProperty(const char* surfaceName,
                                  const char* propertyName,
                                  Int_t np, Double_t* pp, Double_t* values)
{
  Warning("SetMaterialProperty", "Not supported");
}

template <typename K>
Bool_t FastSim<K>::GetTransformation(const TString& volumePath,
                                  TGeoHMatrix& matrix)
{
  Warning("GetTransformation", "Not yet implemented");
  return kFALSE;
}

template <typename K>
Bool_t FastSim<K>::GetShape(const TString& volumePath,
                         TString& shapeType, TArrayD& par)
{
  Warning("GetShape", "Not yet implemented");
  return kFALSE;
}

template <typename K>
Bool_t FastSim<K>::GetMaterial(Int_t imat, TString& name,
                            Double_t& a, Double_t& z, Double_t& density,
                            Double_t& radl, Double_t& inter, TArrayD& par)
{
  Warning("GetMaterial", "Not yet implemented");
  return kFALSE;
}

template <typename K>
Bool_t FastSim<K>::GetMaterial(const TString& volumeName,
                            TString& name, Int_t& imat,
                            Double_t& a, Double_t& z, Double_t& density,
                            Double_t& radl, Double_t& inter, TArrayD& par)
{
  Warning("GetMaterial", "Not yet implemented");
  return kFALSE;
}

template <typename K>
Bool_t FastSim<K>::GetMedium(const TString& volumeName,
                          TString& name, Int_t& imed,
                          Int_t& nmat, Int_t& isvol, Int_t& ifield,
                          Double_t& fieldm, Double_t& tmaxfd, Double_t& stemax,
                          Double_t& deemax, Double_t& epsil, Double_t& stmin,
                          TArrayD& par)
{
  Warning("GetMedium", "Not yet implemented");
  return kFALSE;
}

template <typename K>
void FastSim<K>::WriteEuclid(const char* filnam, const char* topvol,
                          Int_t number, Int_t nlevel)
{
  Warning("WriteEuclid", "Not supported");
}

template <typename K>
void FastSim<K>::SetRootGeometry()
{
  fIsRootGeometrySupported = kTRUE;
}

template <typename K>
void FastSim<K>::SetUserParameters(Bool_t isUserParameters)
{
  Warning("SetUserParameters", "Not supported");
}

template <typename K>
Int_t FastSim<K>::VolId(const char* volName) const
{
  Warning("VolId", "Not yet implemented");
  return -1;
}

template <typename K>
const char* FastSim<K>::VolName(Int_t id) const
{
  Warning("VolName", "Not yet implemented");
  return "";
}

template <typename K>
Int_t FastSim<K>::MediumId(const char* mediumName) const
{
  Warning("MediumId", "Not yet implemented");
  return -1;
}

template <typename K>
Int_t FastSim<K>::NofVolumes() const
{
  Warning("NofVolumes", "Not yet implemented");
  return -1;
}

template <typename K>
Int_t FastSim<K>::VolId2Mate(Int_t id) const
{
  Warning("VolId2Mate", "Not yet implemented");
  return -1;
}

template <typename K>
Int_t FastSim<K>::NofVolDaughters(const char* volName) const
{
  Warning("NofVolDaughters", "Not yet implemented");
  return -1;
}

template <typename K>
const char* FastSim<K>::VolDaughterName(const char* volName, Int_t i) const
{
  Warning("VolDaughterName", "Not yet implemented");
  return "";
}

template <typename K>
Int_t FastSim<K>::VolDaughterCopyNo(const char* volName, Int_t i) const
{
  Warning("VolDaughterCopyNo", "Not yet implemented");
  return -1;
}

template <typename K>
Bool_t FastSim<K>::SetCut(const char* cutName, Double_t cutValue)
{
  Warning("SetCut", "Not yet implemented");
  return kFALSE;
}

template <typename K>
Bool_t FastSim<K>::SetProcess(const char* flagName, Int_t flagValue)
{
  Warning("SetProcess", "Not yet implemented");
  return kFALSE;
}

template <typename K>
Bool_t FastSim<K>::DefineParticle(Int_t pdg, const char* name,
                     TMCParticleType mcType,
                     Double_t mass, Double_t charge, Double_t lifetime)
{
  Warning("DefineParticle", "Not yet implemented");
  return kFALSE;
}

template <typename K>
Bool_t FastSim<K>::DefineParticle(Int_t pdg, const char* name,
                               TMCParticleType mcType,
                               Double_t mass, Double_t charge, Double_t lifetime,
                               const TString& pType, Double_t width,
                               Int_t iSpin, Int_t iParity, Int_t iConjugation,
                               Int_t iIsospin, Int_t iIsospinZ, Int_t gParity,
                               Int_t lepton, Int_t baryon,
                               Bool_t stable, Bool_t shortlived,
                               const TString& subType,
                               Int_t antiEncoding, Double_t magMoment,
                               Double_t excitation)
{
  Warning("DefineParticle", "Not yet implemented");
  return kFALSE;
}

template <typename K>
Bool_t FastSim<K>::DefineIon(const char* name, Int_t Z, Int_t A,
                          Int_t Q, Double_t excEnergy, Double_t mass)
{
  Warning("DefineIon", "Not yet implemented");
  return kFALSE;
}

template <typename K>
Bool_t FastSim<K>::SetDecayMode(Int_t pdg, Float_t bratio[6], Int_t mode[6][3])
{
  Warning("SetDecayMode", "Not yet implemented");
  return kFALSE;
}

template <typename K>
Double_t FastSim<K>::Xsec(char*, Double_t, Int_t, Int_t)
{
  Warning("Xsec", "Not yet implemented");
  return -1.;
}

template <typename K>
Int_t FastSim<K>::IdFromPDG(Int_t pdg) const
{
  Warning("IdFromPDG", "Not yet implemented");
  return -1;
}

template <typename K>
Int_t FastSim<K>::PDGFromId(Int_t id) const
{
  Warning("PDGFromId", "Not yet implemented");
  return -1;
}

template <typename K>
TString FastSim<K>::ParticleName(Int_t pdg) const
{
  Warning("ParticleName", "Not yet implemented");
  return TString();
}

template <typename K>
Double_t FastSim<K>::ParticleMass(Int_t pdg) const
{
  Warning("ParticleMass", "Not yet implemented");
  return -1.;
}

template <typename K>
Double_t FastSim<K>::ParticleCharge(Int_t pdg) const
{
  Warning("ParticleCharge", "Not yet implemented");
  return -1.;
}

template <typename K>
Double_t FastSim<K>::ParticleLifeTime(Int_t pdg) const
{
  Warning("ParticleLifeTime", "Not yet implemented");
  return -1.;
}

template <typename K>
TMCParticleType FastSim<K>::ParticleMCType(Int_t pdg) const
{
  Warning("ParticleMCType", "Not yet implemented");
  return TMCParticleType();
}

template <typename K>
void FastSim<K>::SetMaxStep(Double_t)
{
  Warning("SetMaxStep", "Not yet implemented");
}

template <typename K>
void FastSim<K>::SetMaxNStep(Int_t)
{
  Warning("SetMaxNStep", "Not yet implemented");
}

template <typename K>
void FastSim<K>::SetUserDecay(Int_t pdg)
{
  Warning("SetUserDecay", "Not yet implemented");
}

template <typename K>
void FastSim<K>::ForceDecayTime(Float_t)
{
  Warning("ForceDecayTime", "Not yet implemented");
}

template <typename K>
Int_t FastSim<K>::CurrentVolID(Int_t& copyNo) const
{
  Warning("CurrentVolID", "Not yet implemented");
  return -1;
}

template <typename K>
Int_t FastSim<K>::CurrentVolOffID(Int_t off, Int_t& copyNo) const
{
  Warning("CurrentVolOffID", "Not yet implemented");
  return -1;
}

template <typename K>
const char* FastSim<K>::CurrentVolName() const
{
  Warning("CurrentVolName", "Not yet implemented");
  return "";
}

template <typename K>
const char* FastSim<K>::CurrentVolOffName(Int_t off) const
{
  Warning("CurrentVolOffName", "Not yet implemented");
  return "";
}

template <typename K>
const char* FastSim<K>::CurrentVolPath()
{
  Warning("CurrentVolPath", "Not yet implemented");
  return "";
}

template <typename K>
Bool_t FastSim<K>::CurrentBoundaryNormal(
                    Double_t &x, Double_t &y, Double_t &z) const
{
  Warning("CurrentBoundaryNormal", "Not yet implemented");
  return kFALSE;
}

template <typename K>
Int_t FastSim<K>::CurrentMaterial(Float_t &a, Float_t &z,
                    Float_t &dens, Float_t &radl, Float_t &absl) const
{
  Warning("CurrentMaterial", "Not yet implemented");
  return -1;
}

template <typename K>
Int_t FastSim<K>::CurrentMedium() const
{
  Warning("CurrentMedium", "Not yet implemented");
  return -1;
}

template <typename K>
Int_t FastSim<K>::CurrentEvent() const
{
  Warning("CurrentEvent", "Not yet implemented");
  return -1;
}

template <typename K>
void FastSim<K>::Gmtod(Float_t* xm, Float_t* xd, Int_t iflag)
{
  Warning("Gmtod", "Not yet implemented");
}

template <typename K>
void FastSim<K>::Gmtod(Double_t* xm, Double_t* xd, Int_t iflag)
{
  Warning("Gmtod", "Not yet implemented");
}

template <typename K>
void FastSim<K>::Gdtom(Float_t* xd, Float_t* xm, Int_t iflag)
{
  Warning("Gdtom", "Not yet implemented");
}

template <typename K>
void FastSim<K>::Gdtom(Double_t* xd, Double_t* xm, Int_t iflag)
{
  Warning("Gdtom", "Not yet implemented");
}

template <typename K>
Double_t FastSim<K>::MaxStep() const
{
  Warning("MaxStep", "Not yet implemented");
  return -1.;
}

template <typename K>
Int_t FastSim<K>::GetMaxNStep() const
{
  Warning("GetMaxNStep", "Not yet implemented");
  return -1;
}

template <typename K>
void FastSim<K>::TrackPosition(TLorentzVector& position) const
{
  if(fKernelMode != EKernelMode::kSTEPPING) {
    Warning("TrackPosition", "Kernel not in stepping mode");
    return;
  }
  position = fCurrentParticleStatus->fCurrentPosition;
}

template <typename K>
void FastSim<K>::TrackPosition(Double_t &x, Double_t &y, Double_t &z, Double_t &t) const
{
  if(fKernelMode != EKernelMode::kSTEPPING) {
    Warning("TrackPosition", "Kernel not in stepping mode");
    return;
  }
  ExtractTVector(fCurrentParticleStatus->fCurrentPosition, x, y, z, t);
}

template <typename K>
void FastSim<K>::TrackPosition(Double_t &x, Double_t &y, Double_t &z) const
{
  if(fKernelMode != EKernelMode::kSTEPPING) {
    Warning("TrackPosition", "Kernel not in stepping mode");
    return;
  }
  ExtractTVector(fCurrentParticleStatus->fCurrentPosition, x, y, z);
}

template <typename K>
void FastSim<K>::TrackPosition(Float_t &x, Float_t &y, Float_t &z, Float_t &t) const
{
  if(fKernelMode != EKernelMode::kSTEPPING) {
    Warning("TrackPosition", "Kernel not in stepping mode");
    return;
  }
  ExtractTVector(fCurrentParticleStatus->fCurrentPosition, x, y, z, t);
}

template <typename K>
void FastSim<K>::TrackPosition(Float_t &x, Float_t &y, Float_t &z) const
{
  if(fKernelMode != EKernelMode::kSTEPPING) {
    Warning("TrackPosition", "Kernel not in stepping mode");
    return;
  }
  ExtractTVector(fCurrentParticleStatus->fCurrentPosition, x, y, z);
}

template <typename K>
void FastSim<K>::TrackMomentum(TLorentzVector& momentum) const
{
  if(fKernelMode != EKernelMode::kSTEPPING) {
    Warning("TrackMomentum", "Kernel not in stepping mode");
    return;
  }
  momentum = fCurrentParticleStatus->fCurrentMomentum;
}

template <typename K>
void FastSim<K>::TrackMomentum(Double_t &px, Double_t &py, Double_t &pz, Double_t &etot) const
{
  if(fKernelMode != EKernelMode::kSTEPPING) {
    Warning("TrackMomentum", "Kernel not in stepping mode");
    return;
  }
  ExtractTVector(fCurrentParticleStatus->fCurrentMomentum, px, py, pz, etot);
}

template <typename K>
void FastSim<K>::TrackMomentum(Float_t &px, Float_t &py, Float_t &pz, Float_t &etot) const
{
  if(fKernelMode != EKernelMode::kSTEPPING) {
    Warning("TrackMomentum", "Kernel not in stepping mode");
    return;
  }
  ExtractTVector(fCurrentParticleStatus->fCurrentMomentum, px, py, pz, etot);
}

template <typename K>
Double_t FastSim<K>::TrackStep() const
{
  Warning("TrackStep", "Not yet implemented");
  return -1.;
}

template <typename K>
Double_t FastSim<K>::TrackLength() const
{
  Warning("TrackLength", "Not yet implemented");
  return -1.;
}

template <typename K>
Double_t FastSim<K>::TrackTime() const
{
  Warning("TrackTime", "Not yet implemented");
  return -1.;
}

template <typename K>
Double_t FastSim<K>::Edep() const
{
  Warning("Edep", "Not yet implemented");
  return -1.;
}

template <typename K>
Double_t FastSim<K>::NIELEdep() const
{
  Warning("NIELEdep", "Not yet implemented");
  return -1.;
}

template <typename K>
Int_t FastSim<K>::StepNumber() const
{
  Warning("StepNumber", "Not yet implemented");
  return -1;
}


template <typename K>
Double_t FastSim<K>::TrackWeight() const
{
  Warning("TrackWeight", "Not yet implemented");
  return -1.;
}

template <typename K>
void FastSim<K>::TrackPolarization(Double_t &polX, Double_t &polY, Double_t &polZ) const
{
  Warning("TrackPolarization", "Not yet implemented");
}

template <typename K>
void FastSim<K>::TrackPolarization(TVector3& pol) const
{
  Warning("TrackPolarization", "Not yet implemented");
}

template <typename K>
void FastSim<K>::TrackStatus(TMCParticleStatus& TrackStatus) const
{
  Warning("TrackStatus", "Not yet implemented");
}

template <typename K>
Int_t FastSim<K>::TrackPid() const
{
  Warning("TrackPid", "Not yet implemented");
  return -1;
}

template <typename K>
Double_t FastSim<K>::TrackCharge() const
{
  Warning("TrackCharge", "Not yet implemented");
  return -1.;
}


template <typename K>
Double_t FastSim<K>::TrackMass() const
{
  Warning("TrackMass", "Not yet implemented");
  return -1.;
}


template <typename K>
Double_t FastSim<K>::Etot() const
{
  Warning("Etot", "Not yet implemented");
  return -1.;
}


template <typename K>
Bool_t FastSim<K>::IsNewTrack() const
{
  if(fCurrentParticleStatus->fCurrentStepNumber == 0) {
    return kTRUE;
  }
  return kFALSE;
}

template <typename K>
Bool_t FastSim<K>::IsTrackInside() const
{
  Warning("IsTrackInside", "Not yet implemented");
  return kFALSE;
}

template <typename K>
Bool_t FastSim<K>::IsTrackEntering() const
{
  Warning("IsTrackEntering", "Not yet implemented");
  return kFALSE;
}

template <typename K>
Bool_t FastSim<K>::IsTrackExiting() const
{
  Warning("IsTrackExiting", "Not yet implemented");
  return kFALSE;
}

template <typename K>
Bool_t FastSim<K>::IsTrackOut() const
{
  Warning("IsTrackOut", "Not yet implemented");
  return kFALSE;
}

template <typename K>
Bool_t FastSim<K>::IsTrackDisappeared() const
{
  Warning("IsTrackDisappeared", "Not yet implemented");
  return kFALSE;
}

template <typename K>
Bool_t FastSim<K>::IsTrackStop() const
{
  return fIsTrackStopped;
}

template <typename K>
Bool_t FastSim<K>::IsTrackAlive() const
{
  return !fIsTrackStopped;
}

template <typename K>
Int_t FastSim<K>::NSecondaries() const
{
  Warning("NSecondaries", "Not yet implemented");
  return -1;
}

template <typename K>
void FastSim<K>::GetSecondary(Int_t isec, Int_t& particleId,
                           TLorentzVector& position, TLorentzVector& momentum)
{
  Warning("GetSecondary", "Not yet implemented");
}

template <typename K>
TMCProcess FastSim<K>::ProdProcess(Int_t isec) const
{
  Warning("ProdProcess", "Not yet implemented");
  return TMCProcess();
}

template <typename K>
Int_t FastSim<K>::StepProcesses(TArrayI &proc) const
{
  Warning("StepProcesses", "Not yet implemented");
  return -1;
}

template <typename K>
Bool_t FastSim<K>::SecondariesAreOrdered() const
{
  Warning("SecondariesAreOrdered", "Not yet implemented");
  return kFALSE;
}

template <typename K>
void FastSim<K>::Init()
{
  Warning("Init", "Not yet implemented");
}

template <typename K>
void FastSim<K>::BuildPhysics()
{
  Warning("BuildPhysics", "Not yet implemented");
}

template <typename K>
void FastSim<K>::ProcessEvent(Int_t eventId)
{
  // Generate primaries/get particles from stack
  if(!UseExternalParticleGeneration()) {
    fApplication->GeneratePrimaries();
  }

  fApplication->BeginEvent();

  Int_t trackId = -1;

  while((fCurrentParticle = fMCStack->PopNextTrack(trackId))) {
    if(fMCManagerStack) {
      fCurrentParticleStatus = const_cast<TMCParticleStatus*>(fMCManagerStack->GetParticleStatus(trackId));
    } else {
      fUseTemporaryParticleStatus = kTRUE;
      fCurrentParticleStatus = new TMCParticleStatus();
      fCurrentParticleStatus->InitFromParticle(fCurrentParticle);
      fCurrentParticleStatus->fId = trackId;
    }
    ProcessTrack();
  }
  fApplication->FinishEvent();
}

template <typename K>
void FastSim<K>::ProcessEvent()
{
  Warning("ProcessEvent", "Not yet implemented");
}

template <typename K>
Bool_t FastSim<K>::ProcessRun(Int_t nevent)
{
  if(nevent <= 0) {
    Fatal("ProcessRun", "Number of events must be greater than 0");
  }

  fMCStack = GetStack();
  fMCManagerStack = GetManagerStack();


  for(Int_t i = 0; i < nevent; i++) {
    ProcessEvent(i);
  }
  return kTRUE;
}

template <typename K>
void FastSim<K>::TerminateRun()
{
  Warning("TerminateRun", "Nothing to be done yet");
}

template <typename K>
void FastSim<K>::InitLego()
{
  Warning("InitLego", "Not yet implemented");
}

template <typename K>
void FastSim<K>::SetCollectTracks(Bool_t collectTracks)
{
  Warning("SetCollectTracks", "Not yet implemented");
}

template <typename K>
Bool_t FastSim<K>::IsCollectTracks() const
{
  Warning("IsCollectTracks", "Not yet implemented");
  return kFALSE;
}


//ClassImp(vmcfastsim::base::FastSim<K>)
