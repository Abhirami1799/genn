#include "currentSource.h"

// Standard includes
#include <algorithm>
#include <cmath>

// GeNN includes
#include "gennUtils.h"

//------------------------------------------------------------------------
// CurrentSource
//------------------------------------------------------------------------
void CurrentSource::setVarLocation(const std::string &varName, VarLocation loc)
{
    m_VarLocation[getCurrentSourceModel()->getVarIndex(varName)] = loc;
}
//----------------------------------------------------------------------------
void CurrentSource::setVarImplementation(const std::string &varName, VarImplementation impl)
{
    m_VarImplementation[getCurrentSourceModel()->getVarIndex(varName)] = impl;
}
//----------------------------------------------------------------------------
void CurrentSource::setExtraGlobalParamLocation(const std::string &paramName, VarLocation loc)
{
    const size_t extraGlobalParamIndex = getCurrentSourceModel()->getExtraGlobalParamIndex(paramName);
    if(!Utils::isTypePointer(getCurrentSourceModel()->getExtraGlobalParams()[extraGlobalParamIndex].type)) {
        throw std::runtime_error("Only extra global parameters with a pointer type have a location");
    }
    m_ExtraGlobalParamLocation[extraGlobalParamIndex] = loc;
}
//----------------------------------------------------------------------------
VarLocation CurrentSource::getVarLocation(const std::string &varName) const
{
    return m_VarLocation[getCurrentSourceModel()->getVarIndex(varName)];
}
//----------------------------------------------------------------------------
VarImplementation CurrentSource::getVarImplementation(const std::string &varName) const
{
    return m_VarImplementation[getCurrentSourceModel()->getVarIndex(varName)];
}
//----------------------------------------------------------------------------
VarLocation CurrentSource::getExtraGlobalParamLocation(const std::string &varName) const
{
    return m_ExtraGlobalParamLocation[getCurrentSourceModel()->getExtraGlobalParamIndex(varName)];
}
//----------------------------------------------------------------------------
CurrentSource::CurrentSource(const std::string &name, const CurrentSourceModels::Base *currentSourceModel,
                             const std::vector<double> &params, const std::vector<Models::VarInit> &varInitialisers,
                             VarLocation defaultVarLocation, VarLocation defaultExtraGlobalParamLocation)
:   m_Name(name), m_CurrentSourceModel(currentSourceModel), m_VarLocation(params.size() + varInitialisers.size(), defaultVarLocation),
    m_ExtraGlobalParamLocation(currentSourceModel->getExtraGlobalParams().size(), defaultExtraGlobalParamLocation)
{
    // Populate combined variable initialisers and implementations from legacy parameters and initialisers
    Utils::initialiseLegacyImplementation(params, varInitialisers, m_VarInitialisers, m_VarImplementation);
}
//----------------------------------------------------------------------------
CurrentSource::CurrentSource(const std::string &name, const CurrentSourceModels::Base *currentSourceModel,
                             const std::vector<Models::VarInit> &varInitialisers,
                             VarLocation defaultVarLocation, VarLocation defaultExtraGlobalParamLocation)
:   m_Name(name), m_CurrentSourceModel(currentSourceModel), m_VarInitialisers(varInitialisers),
    m_VarLocation(varInitialisers.size(), defaultVarLocation),
    m_ExtraGlobalParamLocation(currentSourceModel->getExtraGlobalParams().size(), defaultExtraGlobalParamLocation)
{
    if(!currentSourceModel->getParamNames().empty()) {
        throw std::runtime_error("Populations using current source models with legacy 'parameters' cannot be added with this functions");
    }

    // Automatically determine default implementation for variables
    Utils::autoDetermineImplementation(varInitialisers, currentSourceModel->getCombinedVars(), m_VarImplementation);
}
//----------------------------------------------------------------------------
void CurrentSource::initDerivedParams(double dt)
{
    // Calculate derived parameter values
    Utils::calcDerivedParamVal(getCurrentSourceModel(), m_VarInitialisers, m_VarImplementation, dt, m_DerivedParams);

    // Initialise derived parameters for variable initialisers
    for(auto &v : m_VarInitialisers) {
        v.initDerivedParams(dt);
    }
}
//----------------------------------------------------------------------------
bool CurrentSource::isSimRNGRequired() const
{
    // Returns true if any parts of the current source code require an RNG
    if(Utils::isRNGRequired(getCurrentSourceModel()->getInjectionCode())) {
        return true;
    }

    return false;
}
//----------------------------------------------------------------------------
bool CurrentSource::isInitRNGRequired() const
{
    // If initialising the neuron variables require an RNG, return true
    if(Utils::isInitRNGRequired(getVarInitialisers())) {
        return true;
    }

    return false;
}
