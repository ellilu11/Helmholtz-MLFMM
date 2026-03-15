#pragma once

#include "config.h"
#include "math.h"

class Source;

using SrcVec = std::vector<std::shared_ptr<Source>>;

class Source {

public:
    Source() = default;
    
    Source(size_t iSrc) : iSrc(iSrc), rval(0.0) {};

    size_t getIdx() const { return iSrc; }

    void setIdx(size_t iSrc) { this->iSrc = iSrc; }

    void addToRval(cmplx val) { rval += val; }

    template <typename T>
    bool isSrcType() const { return typeid(*this) == typeid(T); }

    virtual vec3d getCenter() const = 0;

    virtual cmplx getVoltage() = 0;

    virtual std::pair<vec3cd,vec3cd> getRadsAlongDir(const vec3d&, const vec3d&) const = 0;

    virtual vec3cd getFarAlongDir(const vec3d&) const = 0;

    virtual cmplx getIntegratedEFIE(const std::shared_ptr<Source>) const = 0;

    virtual cmplx getIntegratedMFIE(const std::shared_ptr<Source>) const = 0;

    virtual double getIntegratedMass(const std::shared_ptr<Source>) const = 0;

protected:
    size_t iSrc;
    cmplx rval;
};