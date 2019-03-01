// Some license stuff


#ifndef FAST_SIM_KERNEL_H
#define FAST_SIM_KERNEL_H

namespace vmcfastsim {

namespace base {

class FastSimTrack;

class VFastSimKernel
{
  public:
    VFastSimKernel()
      : mTrack(nullptr), mIsValid(false)
    {}
    virtual ~VFastSimKernel() = default;

    void SetTrack()
    {
      // Nothing to be done yet
    }

    void SetValidCheck()
    {
      // Nothing to be done yet
    }

    bool IsValid() const
    {
      return mIsValid;
    }

  protected:
    // Can be called to tell that a step was made
    void Stepping()
    {
      // Nothing to be done yet
    }

    /// Use to forward newly generated particle
    void Push(FastSimTrack* track)
    {
      // Nothing to be done yet
    }

  protected:
    FastSimTrack* mTrack;

  private:
    bool mIsValid;

};


// Template parameter I defining an interface calling user hooks and communicate
// with the user stack
template <class U>
class FastSimKernelImpl : public VFastSimKernel
{
  public:
    FastSimKernelImpl()
      : VFastSimKernel()
    {}
    virtual ~FastSimKernelImpl() = default;


    bool Process() override final
    {
      IsValid();
      bool returnValue = static_cast<U*>(this)->ProcessImpl();
      //Post();
      return returnValue;
    }

    // This initializes the kernel and initializes also the interface if not yet
    // done
    void Initialize() override final
    {
      static_cast<U*>(this)->InitializeImpl();
    }

    void Stop() override final
    {
      static_cast<U*>(this)->StopImpl();
    }

    virtual bool ProcessImpl() = 0;
    virtual void InitializeImpl() = 0;
    virtual void StopImpl() = 0;

};

} // namespace base

} // namespace vmcfastsim



#endif /* FAST_SIM_KERNEL_H */
