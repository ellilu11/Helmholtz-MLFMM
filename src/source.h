#pragma once

#include <filesystem>
#include "config.h"
#include "excite.h"
#include "math.h"
#include "phys.h"

class Source;

using SrcVec = std::vector<std::shared_ptr<Source>>;

class Source {

public:
    Source() = default;
    
    Source(std::shared_ptr<Exc::PlaneWave> Einc, size_t iSrc)
        : Einc(std::move(Einc)), iSrc(iSrc), voltage(0.0), rval(0.0)
    {};

    cmplx getVoltage() const { return voltage; }

    size_t getIdx() const { return iSrc; }

    void setIdx(size_t iSrc) { this->iSrc = iSrc; }

    void addToRval(cmplx val) { rval += val; }

    template <typename T>
    bool isSrcType() const { return typeid(*this) == typeid(T); }

    virtual vec3d getCenter() const = 0;

    virtual void buildVoltage() = 0;

    virtual vec3cd getRadAlongDir(const vec3d&, const vec3d&) const = 0;

    virtual vec3cd getFarAlongDir(const vec3d&) const = 0;

    virtual cmplx getIntegratedEFIE(const std::shared_ptr<Source>) const = 0;

    virtual cmplx getIntegratedMFIE(const std::shared_ptr<Source>) const = 0;

protected:
    std::shared_ptr<Exc::PlaneWave> Einc;
    cmplx voltage;
    size_t iSrc;

    cmplx rval;
};