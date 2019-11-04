# Fast sim base class based on `TVirtualMC`

**Disclaimer:** This is not a class yet fully optimised to serve as a base for state-of-the-art fast simulation implementation. However, it already implements all pure virtual methods of `TVirtualMC` and makes it hence easy, to quickly write a custom VMC. For an example see [VMCFastShower](https://github.com/benedikt-voelkel/FastShower).

## Installation
The small package is set up using CMake in a dedicated build directory via

```bash
mkdir $BUILD_DIR
cd $BUILD_DIR
cmake $SOURCE_DIR -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR -DROOT_DIR=$ROOT_CMAKE_CONFIG -DVMC_DIR=$VMC_INSTALL_DIR
make && make install
```

## Use this base class

A custom user implementation using this as a base looks like the following
```cpp
class FastShower : public::vmcfastsim::base::FastSim<FastShower>
{
    FastShower(...);
    // ...
    bool Process() final
    {
        // ...
    }

    void Stop() final
    {
        // Nothing to do yet
    }

    // ...
};
```

This is bsaically taken from [VMCFastShower](https://github.com/benedikt-voelkel/FastShower/blob/master/include/FastShower.h)
