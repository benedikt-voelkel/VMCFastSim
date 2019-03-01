// Some license stuff


#ifndef FAST_SIM_GEOMETRY_INTERFACE_H
#define FAST_SIM_GEOMETRY_INTERFACE_H

namespace vmcfastsim {

namespace geometry {



class FastSimGeometryInterfaceBase
{

  FastSimGeometryInterfaceBase() = default;
  virtual ~FastSimGeometryInterfaceBase() = default;


  protected:
    void InitializeTrack(TMCParticleStatus* particleStatus)
    {
      fCurrentParticleStatus = particleStatus;
    }

  protected:
    TMCParticleStatus* fCurrentParticleStatus;
}


class FastSimGeometryInterfaceROOT : public FastSimGeometryInterfaceBase
{
  public:
    FastSimGeometryInterfaceROOT()
      : FastSimGeometryInterfaceBase()
    {}

    virtual ~FastSimGeometryInterfaceROOT() = default;


    /// Update according to fCurrentParticleStatus
    void Update()
    {

    }





};


} // namespace geometry

} // namespace vmcfastsim


#endif /* FAST_SIM_GEOMETRY_INTERFACE_H */
