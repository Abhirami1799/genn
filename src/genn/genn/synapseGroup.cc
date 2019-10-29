#include "synapseGroup.h"

// Standard includes
#include <algorithm>
#include <cmath>
#include <iostream>

// GeNN includes
#include "neuronGroupInternal.h"
#include "gennUtils.h"

// ------------------------------------------------------------------------
// SynapseGroup
// ------------------------------------------------------------------------
void SynapseGroup::setWUVarLocation(const std::string &varName, VarLocation loc)
{
    m_WUVarLocation[getWUModel()->getVarIndex(varName)] = loc;
}
//----------------------------------------------------------------------------
void SynapseGroup::setWUVarImplementation(const std::string &varName, VarImplementation impl)
{
    m_WUVarImplementation[getWUModel()->getVarIndex(varName)] = impl;
}
//----------------------------------------------------------------------------
void SynapseGroup::setWUPreVarLocation(const std::string &varName, VarLocation loc)
{
    m_WUPreVarLocation[getWUModel()->getPreVarIndex(varName)] = loc;
}
//----------------------------------------------------------------------------
void SynapseGroup::setWUPreVarImplementation(const std::string &varName, VarImplementation impl)
{
    m_WUPreVarImplementation[getWUModel()->getPreVarIndex(varName)] = impl;
}
//----------------------------------------------------------------------------
void SynapseGroup::setWUPostVarLocation(const std::string &varName, VarLocation loc)
{
    m_WUPostVarLocation[getWUModel()->getPostVarIndex(varName)] = loc;
}
//----------------------------------------------------------------------------
void SynapseGroup::setWUPostVarImplementation(const std::string &varName, VarImplementation impl)
{
    m_WUPostVarImplementation[getWUModel()->getPostVarIndex(varName)] = impl;
}
//----------------------------------------------------------------------------
void SynapseGroup::setWUExtraGlobalParamLocation(const std::string &paramName, VarLocation loc)
{
    const size_t extraGlobalParamIndex = getWUModel()->getExtraGlobalParamIndex(paramName);
    if(!Utils::isTypePointer(getWUModel()->getExtraGlobalParams()[extraGlobalParamIndex].type)) {
        throw std::runtime_error("Only extra global parameters with a pointer type have a location");
    }
    m_WUExtraGlobalParamLocation[extraGlobalParamIndex] = loc;
}
//----------------------------------------------------------------------------
void SynapseGroup::setPSVarLocation(const std::string &varName, VarLocation loc)
{
    m_PSVarLocation[getPSModel()->getVarIndex(varName)] = loc;
}
//----------------------------------------------------------------------------
void SynapseGroup::setPSVarImplementation(const std::string &varName, VarImplementation impl)
{
    m_PSVarImplementation[getPSModel()->getVarIndex(varName)] = impl;
}
//----------------------------------------------------------------------------
void SynapseGroup::setPSExtraGlobalParamLocation(const std::string &paramName, VarLocation loc)
{
    const size_t extraGlobalParamIndex = getPSModel()->getExtraGlobalParamIndex(paramName);
    if(!Utils::isTypePointer(getPSModel()->getExtraGlobalParams()[extraGlobalParamIndex].type)) {
        throw std::runtime_error("Only extra global parameters with a pointer type have a location");
    }
    m_PSExtraGlobalParamLocation[extraGlobalParamIndex] = loc;
}
//----------------------------------------------------------------------------
void SynapseGroup::setSparseConnectivityExtraGlobalParamLocation(const std::string &paramName, VarLocation loc)
{
    const size_t extraGlobalParamIndex = m_ConnectivityInitialiser.getSnippet()->getExtraGlobalParamIndex(paramName);
    if(!Utils::isTypePointer(m_ConnectivityInitialiser.getSnippet()->getExtraGlobalParams()[extraGlobalParamIndex].type)) {
        throw std::runtime_error("Only extra global parameters with a pointer type have a location");
    }
    m_ConnectivityExtraGlobalParamLocation[extraGlobalParamIndex] = loc;
}
//----------------------------------------------------------------------------
void SynapseGroup::setMaxConnections(unsigned int maxConnections)
{
    if(m_MatrixConnectivity == SynapseMatrixConnectivity::SPARSE) {
        if(m_ConnectivityInitialiser.getSnippet()->getCalcMaxRowLengthFunc()) {
            throw std::runtime_error("setMaxConnections: Synapse group already has max connections defined by connectivity initialisation snippet.");
        }
        
        m_MaxConnections = maxConnections;
    }
    else {
        throw std::runtime_error("setMaxConnections: Synapse group is densely connected. Setting max connections is not required in this case.");
    }
}
//----------------------------------------------------------------------------
void SynapseGroup::setMaxSourceConnections(unsigned int maxConnections)
{
    if(m_MatrixConnectivity == SynapseMatrixConnectivity::SPARSE) {
        if(m_ConnectivityInitialiser.getSnippet()->getCalcMaxColLengthFunc()) {
            throw std::runtime_error("setMaxSourceConnections: Synapse group already has max source connections defined by connectivity initialisation snippet.");
        }

        m_MaxSourceConnections = maxConnections;
    }
    else {
        throw std::runtime_error("setMaxSourceConnections: Synapse group is densely connected. Setting max connections is not required in this case.");
    }
}
//----------------------------------------------------------------------------
void SynapseGroup::setMaxDendriticDelayTimesteps(unsigned int maxDendriticDelayTimesteps)
{
    // **TODO** constraints on this
    m_MaxDendriticDelayTimesteps = maxDendriticDelayTimesteps;
}
//----------------------------------------------------------------------------
void SynapseGroup::setSpanType(SpanType spanType)
{
    if(m_MatrixConnectivity == SynapseMatrixConnectivity::SPARSE || m_MatrixConnectivity == SynapseMatrixConnectivity::PROCEDURAL) {
        m_SpanType = spanType;
    }
    else {
        throw std::runtime_error("setSpanType: This function can only be used on synapse groups with sparse connectivity.");
    }
}
//----------------------------------------------------------------------------
void SynapseGroup::setNumThreadsPerSpike(unsigned int numThreadsPerSpike)
{
    if(m_SpanType == SpanType::PRESYNAPTIC) {
        m_NumThreadsPerSpike = numThreadsPerSpike;
    }
    else {
        throw std::runtime_error("setSpanType: This function can only be used on synapse groups with a presynaptic span type.");
    }
}
//----------------------------------------------------------------------------
void SynapseGroup::setBackPropDelaySteps(unsigned int timesteps)
{
    m_BackPropDelaySteps = timesteps;

    m_TrgNeuronGroup->checkNumDelaySlots(m_BackPropDelaySteps);
}
//----------------------------------------------------------------------------
void SynapseGroup::setNarrowSparseIndEnabled(bool enabled)
{
    if (getMatrixConnectivity() == SynapseMatrixConnectivity::SPARSE) {
        m_NarrowSparseIndEnabled = enabled;
    }
    else {
        throw std::runtime_error("setNarrowSparseIndEnabled: This function can only be used on synapse groups with sparse connectivity.");
    }
}
//----------------------------------------------------------------------------
int SynapseGroup::getClusterHostID() const
{
    return m_TrgNeuronGroup->getClusterHostID();
}
//----------------------------------------------------------------------------
bool SynapseGroup::isTrueSpikeRequired() const
{
    return !getWUModel()->getSimCode().empty();
}
//----------------------------------------------------------------------------
bool SynapseGroup::isSpikeEventRequired() const
{
     return !getWUModel()->getEventCode().empty();
}
//----------------------------------------------------------------------------
bool SynapseGroup::isZeroCopyEnabled() const
{
    // If there are any postsynaptic variables implemented in zero-copy mode return true
    if(std::any_of(m_PSVarLocation.begin(), m_PSVarLocation.end(),
        [](VarLocation loc){ return (loc & VarLocation::ZERO_COPY); }))
    {
        return true;
    }

    // If there are any weight update variables implemented in zero-copy mode return true
    if(std::any_of(m_WUVarLocation.begin(), m_WUVarLocation.end(),
        [](VarLocation loc){ return (loc & VarLocation::ZERO_COPY); }))
    {
        return true;
    }

    // If there are any weight update variables implemented in zero-copy mode return true
    if(std::any_of(m_WUPreVarLocation.begin(), m_WUPreVarLocation.end(),
        [](VarLocation loc){ return (loc & VarLocation::ZERO_COPY); }))
    {
        return true;
    }

    // If there are any weight update variables implemented in zero-copy mode return true
    if(std::any_of(m_WUPostVarLocation.begin(), m_WUPostVarLocation.end(),
        [](VarLocation loc){ return (loc & VarLocation::ZERO_COPY); }))
    {
        return true;
    }

    return false;
}
//----------------------------------------------------------------------------
VarLocation SynapseGroup::getWUVarLocation(const std::string &var) const
{
    return m_WUVarLocation[getWUModel()->getVarIndex(var)];
}
//----------------------------------------------------------------------------
VarLocation SynapseGroup::getWUPreVarLocation(const std::string &var) const
{
    return m_WUPreVarLocation[getWUModel()->getPreVarIndex(var)];
}
//----------------------------------------------------------------------------
VarLocation SynapseGroup::getWUPostVarLocation(const std::string &var) const
{
    return m_WUPostVarLocation[getWUModel()->getPostVarIndex(var)];
}
//----------------------------------------------------------------------------
VarImplementation SynapseGroup::getWUVarImplementation(const std::string &var) const
{
    return m_WUVarImplementation[getWUModel()->getVarIndex(var)];
}
//----------------------------------------------------------------------------
VarImplementation SynapseGroup::getWUPreVarImplementation(const std::string &var) const
{
    return m_WUPreVarImplementation[getWUModel()->getPreVarIndex(var)];
}
//----------------------------------------------------------------------------
VarImplementation SynapseGroup::getWUPostVarImplementation(const std::string &var) const
{
    return m_WUPostVarImplementation[getWUModel()->getPostVarIndex(var)];
}
//----------------------------------------------------------------------------
VarLocation SynapseGroup::getWUExtraGlobalParamLocation(const std::string &paramName) const
{
    return m_WUExtraGlobalParamLocation[getWUModel()->getExtraGlobalParamIndex(paramName)];
}
//----------------------------------------------------------------------------
VarLocation SynapseGroup::getPSVarLocation(const std::string &var) const
{
    return m_PSVarLocation[getPSModel()->getVarIndex(var)];
}
//----------------------------------------------------------------------------
VarImplementation SynapseGroup::getPSVarImplementation(const std::string &var) const
{
    return m_PSVarImplementation[getPSModel()->getVarIndex(var)];
}
//----------------------------------------------------------------------------
VarLocation SynapseGroup::getPSExtraGlobalParamLocation(const std::string &paramName) const
{
    return m_PSExtraGlobalParamLocation[getPSModel()->getExtraGlobalParamIndex(paramName)];
}
//----------------------------------------------------------------------------
VarLocation SynapseGroup::getSparseConnectivityExtraGlobalParamLocation(const std::string &paramName) const
{
    return m_ConnectivityExtraGlobalParamLocation[m_ConnectivityInitialiser.getSnippet()->getExtraGlobalParamIndex(paramName)];
}
//----------------------------------------------------------------------------
bool SynapseGroup::isDendriticDelayRequired() const
{
    // If addToInSynDelay function is used in sim code, return true
    if(getWUModel()->getSimCode().find("$(addToInSynDelay") != std::string::npos) {
        return true;
    }

    // If addToInSynDelay function is used in synapse dynamics, return true
    if(getWUModel()->getSynapseDynamicsCode().find("$(addToInSynDelay") != std::string::npos) {
        return true;
    }

    return false;
}
//----------------------------------------------------------------------------
bool SynapseGroup::isProceduralConnectivityRNGRequired() const
{
    return ((m_MatrixType & SynapseMatrixConnectivity::PROCEDURAL) &&
            Utils::isRNGRequired(m_ConnectivityInitialiser.getSnippet()->getRowBuildCode()));
}
//----------------------------------------------------------------------------
bool SynapseGroup::isPSInitRNGRequired() const
{
    // If initialising the postsynaptic variables require an RNG, return true
    return Utils::isInitRNGRequired(m_PSVarInitialisers);
}
//----------------------------------------------------------------------------
bool SynapseGroup::isWUInitRNGRequired() const
{
    // If initialising the weight update variables require an RNG, return true
    if(Utils::isInitRNGRequired(m_WUVarInitialisers)) {
        return true;
    }

    // Return true if matrix has sparse or bitmask connectivity and an RNG is required to initialise connectivity
    return (((m_MatrixType & SynapseMatrixConnectivity::SPARSE) || (m_MatrixType & SynapseMatrixConnectivity::BITMASK))
            && Utils::isRNGRequired(m_ConnectivityInitialiser.getSnippet()->getRowBuildCode()));
}
//----------------------------------------------------------------------------
bool SynapseGroup::isWUVarInitRequired() const
{
    // Loop through variable's initialisers and implementation
    auto varInit = m_WUVarInitialisers.cbegin();
    auto varImpl = m_WUVarImplementation.cbegin();
    for(;varInit != m_WUVarInitialisers.cend(); varInit++, varImpl++) {
        // If variable is individually implemented and has some initialisation code, return true
        if(*varImpl == VarImplementation::INDIVIDUAL && !varInit->getSnippet()->getCode().empty()) {
            return true;
        }
    }

    return false;
}
//----------------------------------------------------------------------------
bool SynapseGroup::isSparseConnectivityInitRequired() const
{
    // Return true if the matrix type is sparse or bitmask and there is code to initialise sparse connectivity
    return (((m_MatrixType & SynapseMatrixConnectivity::SPARSE) || (m_MatrixType & SynapseMatrixConnectivity::BITMASK))
            && !getConnectivityInitialiser().getSnippet()->getRowBuildCode().empty());
}
//----------------------------------------------------------------------------
SynapseGroup::SynapseGroup(const std::string name, SynapseMatrixType matrixType, unsigned int delaySteps,
                           const WeightUpdateModels::Base *wu, const std::vector<double> &wuParams, const std::vector<Models::VarInit> &wuVarInitialisers, const std::vector<Models::VarInit> &wuPreVarInitialisers, const std::vector<Models::VarInit> &wuPostVarInitialisers,
                           const PostsynapticModels::Base *ps, const std::vector<double> &psParams, const std::vector<Models::VarInit> &psVarInitialisers,
                           NeuronGroupInternal *srcNeuronGroup, NeuronGroupInternal *trgNeuronGroup,
                           const InitSparseConnectivitySnippet::Init &connectivityInitialiser,
                           VarLocation defaultVarLocation, VarLocation defaultExtraGlobalParamLocation,
                           VarLocation defaultSparseConnectivityLocation, bool defaultNarrowSparseIndEnabled)
:   m_Name(name), m_SpanType(SpanType::POSTSYNAPTIC), m_NumThreadsPerSpike(1), m_DelaySteps(delaySteps), m_BackPropDelaySteps(0),
    m_MaxDendriticDelayTimesteps(1), m_MatrixConnectivity(getSynapseMatrixConnectivity(matrixType)),
    m_SrcNeuronGroup(srcNeuronGroup), m_TrgNeuronGroup(trgNeuronGroup), m_EventThresholdReTestRequired(false),
    m_NarrowSparseIndEnabled(defaultNarrowSparseIndEnabled), m_InSynLocation(defaultVarLocation),  m_DendriticDelayLocation(defaultVarLocation),
    m_WUModel(wu), m_PSModel(ps), m_WUVarLocation(wuParams.size() + wuVarInitialisers.size(), defaultVarLocation),
    m_WUPreVarLocation(wuPreVarInitialisers.size(), defaultVarLocation), m_WUPostVarLocation(wuPostVarInitialisers.size(), defaultVarLocation),
    m_WUExtraGlobalParamLocation(wu->getExtraGlobalParams().size(), defaultExtraGlobalParamLocation),
    m_PSVarLocation(psParams.size() + psVarInitialisers.size(), defaultVarLocation),
    m_PSExtraGlobalParamLocation(ps->getExtraGlobalParams().size(), defaultExtraGlobalParamLocation),
    m_ConnectivityInitialiser(connectivityInitialiser), m_SparseConnectivityLocation(defaultSparseConnectivityLocation),
    m_ConnectivityExtraGlobalParamLocation(connectivityInitialiser.getSnippet()->getExtraGlobalParams().size(), defaultExtraGlobalParamLocation), m_PSModelTargetName(name)
{
    // Populate combined weight update model variable initialisers and implementations from legacy parameters and initialisers
    Utils::initialiseLegacyImplementation(wuParams, wuVarInitialisers, m_WUVarInitialisers, m_WUVarImplementation,
        (matrixType & SynapseMatrixWeight::INDIVIDUAL) ? VarImplementation::INDIVIDUAL : VarImplementation::GLOBAL);
    Utils::initialiseLegacyImplementation({}, wuPreVarInitialisers, m_WUPreVarInitialisers, m_WUPreVarImplementation);
    Utils::initialiseLegacyImplementation({}, wuPostVarInitialisers, m_WUPostVarInitialisers, m_WUPostVarImplementation);

    // Populate combined weight update model variable initialisers and implementations from legacy parameters and initialisers
    // **NOTE** because legacy parameters are accessible everywhere, add them to all three
    Utils::initialiseLegacyImplementation(psParams, psVarInitialisers, m_PSVarInitialisers, m_PSVarImplementation,
        (matrixType & SynapseMatrixWeight::INDIVIDUAL_PSM) ? VarImplementation::INDIVIDUAL : VarImplementation::GLOBAL);

    setDefaultMaxConnections();
}
//----------------------------------------------------------------------------
SynapseGroup::SynapseGroup(const std::string name, SynapseMatrixConnectivity matrixConnectivity, unsigned int delaySteps,
                           const WeightUpdateModels::Base *wu, const std::vector<Models::VarInit> &wuVarInitialisers, const std::vector<Models::VarInit> &wuPreVarInitialisers, const std::vector<Models::VarInit> &wuPostVarInitialisers,
                           const PostsynapticModels::Base *ps, const std::vector<Models::VarInit> &psVarInitialisers,
                           NeuronGroupInternal *srcNeuronGroup, NeuronGroupInternal *trgNeuronGroup,
                           const InitSparseConnectivitySnippet::Init &connectivityInitialiser,
                           VarLocation defaultVarLocation, VarLocation defaultExtraGlobalParamLocation,
                           VarLocation defaultSparseConnectivityLocation, bool defaultNarrowSparseIndEnabled)
:   m_Name(name), m_SpanType(SpanType::POSTSYNAPTIC), m_NumThreadsPerSpike(1), m_DelaySteps(delaySteps), m_BackPropDelaySteps(0),
    m_MaxDendriticDelayTimesteps(1), m_MatrixConnectivity(matrixConnectivity),  m_SrcNeuronGroup(srcNeuronGroup), m_TrgNeuronGroup(trgNeuronGroup),
    m_EventThresholdReTestRequired(false), m_NarrowSparseIndEnabled(defaultNarrowSparseIndEnabled), 
    m_InSynLocation(defaultVarLocation),  m_DendriticDelayLocation(defaultVarLocation),
    m_WUModel(wu), m_WUVarInitialisers(wuVarInitialisers), m_WUPreVarInitialisers(wuPreVarInitialisers),
    m_WUPostVarInitialisers(wuPostVarInitialisers), m_PSModel(ps), m_PSVarInitialisers(psVarInitialisers),
    m_WUVarLocation(wuVarInitialisers.size(), defaultVarLocation), m_WUPreVarLocation(wuPreVarInitialisers.size(), defaultVarLocation),
    m_WUPostVarLocation(wuPostVarInitialisers.size(), defaultVarLocation), m_WUExtraGlobalParamLocation(wu->getExtraGlobalParams().size(), defaultExtraGlobalParamLocation),
    m_PSVarLocation(psVarInitialisers.size(), defaultVarLocation), m_PSExtraGlobalParamLocation(ps->getExtraGlobalParams().size(), defaultExtraGlobalParamLocation),
    m_ConnectivityInitialiser(connectivityInitialiser), m_SparseConnectivityLocation(defaultSparseConnectivityLocation),
    m_ConnectivityExtraGlobalParamLocation(connectivityInitialiser.getSnippet()->getExtraGlobalParams().size(), defaultExtraGlobalParamLocation), m_PSModelTargetName(name)
{
    if(!wu->getParamNames().empty()) {
        throw std::runtime_error("Populations using weight update models with legacy 'parameters' cannot be added with this functions");
    }

    if(!ps->getParamNames().empty()) {
        throw std::runtime_error("Populations using postsynaptic models with legacy 'parameters' cannot be added with this functions");
    }

    // Automatically determine default implementation for weight update model variables
    Utils::autoDetermineImplementation(m_WUVarInitialisers, wu->getCombinedVars(), m_WUVarImplementation);
    Utils::autoDetermineImplementation(m_WUPreVarInitialisers, wu->getPreVars(), m_WUPreVarImplementation);
    Utils::autoDetermineImplementation(m_WUPostVarInitialisers, wu->getPostVars(), m_WUPostVarImplementation);

    // Automatically determine default implementation for postsynaptic model variables
    Utils::autoDetermineImplementation(m_PSVarInitialisers, ps->getCombinedVars(), m_PSVarImplementation);

    setDefaultMaxConnections();
}
//----------------------------------------------------------------------------
void SynapseGroup::initDerivedParams(double dt)
{
    // Calculate derived parameter values
    Utils::calcDerivedParamVal(getWUModel(), m_WUVarInitialisers, m_WUVarImplementation, dt, m_WUDerivedParams);
    Utils::calcDerivedParamVal(getPSModel(), m_PSVarInitialisers, m_PSVarImplementation, dt, m_PSDerivedParams);

    // Initialise derived parameters for WU variable initialisers
    for(auto &v : m_WUVarInitialisers) {
        v.initDerivedParams(dt);
    }

    // Initialise derived parameters for PSM variable initialisers
    for(auto &v : m_PSVarInitialisers) {
        v.initDerivedParams(dt);
    }

    // Initialise derived parameters for WU presynaptic variable initialisers
    for(auto &v : m_WUPreVarInitialisers) {
        v.initDerivedParams(dt);
    }
    
    // Initialise derived parameters for WU postsynaptic variable initialisers
    for(auto &v : m_WUPostVarInitialisers) {
        v.initDerivedParams(dt);
    }

    // Initialise any derived connectivity initialiser parameters
    m_ConnectivityInitialiser.initDerivedParams(dt);
}
//----------------------------------------------------------------------------
std::string SynapseGroup::getPresynapticAxonalDelaySlot(const std::string &devPrefix) const
{
    assert(getSrcNeuronGroup()->isDelayRequired());

    if(getDelaySteps() == 0) {
        return devPrefix + "spkQuePtr" + getSrcNeuronGroup()->getName();
    }
    else {
        return "((" + devPrefix + "spkQuePtr" + getSrcNeuronGroup()->getName() + " + " + std::to_string(getSrcNeuronGroup()->getNumDelaySlots() - getDelaySteps()) + ") % " + std::to_string(getSrcNeuronGroup()->getNumDelaySlots()) + ")";
    }
}
//----------------------------------------------------------------------------
std::string SynapseGroup::getPostsynapticBackPropDelaySlot(const std::string &devPrefix) const
{
    assert(getTrgNeuronGroup()->isDelayRequired());

    if(getBackPropDelaySteps() == 0) {
        return devPrefix + "spkQuePtr" + getTrgNeuronGroup()->getName();
    }
    else {
        return "((" + devPrefix + "spkQuePtr" + getTrgNeuronGroup()->getName() + " + " + std::to_string(getTrgNeuronGroup()->getNumDelaySlots() - getBackPropDelaySteps()) + ") % " + std::to_string(getTrgNeuronGroup()->getNumDelaySlots()) + ")";
    }
}
//----------------------------------------------------------------------------
std::string SynapseGroup::getDendriticDelayOffset(const std::string &devPrefix, const std::string &offset) const
{
    assert(isDendriticDelayRequired());

    if(offset.empty()) {
        return "(" + devPrefix + "denDelayPtr" + getPSModelTargetName() + " * " + std::to_string(getTrgNeuronGroup()->getNumNeurons()) + ") + ";
    }
    else {
        return "(((" + devPrefix + "denDelayPtr" + getPSModelTargetName() + " + " + offset + ") % " + std::to_string(getMaxDendriticDelayTimesteps()) + ") * " + std::to_string(getTrgNeuronGroup()->getNumNeurons()) + ") + ";
    }
}
//----------------------------------------------------------------------------
std::string SynapseGroup::getSparseIndType() const
{
    // If narrow sparse inds are enabled
    if(m_NarrowSparseIndEnabled) {
        // If number of target neurons can be represented using a uint8, use this type
        const unsigned int numTrgNeurons = getTrgNeuronGroup()->getNumNeurons();
        if(numTrgNeurons <= std::numeric_limits<uint8_t>::max()) {
            return "uint8_t";
        }
        // Otherwise, if they can be represented as a uint16, use this type
        else if(numTrgNeurons <= std::numeric_limits<uint16_t>::max()) {
            return "uint16_t";
        }
    }

    // Otherwise, use 32-bit int
    return "uint32_t";

}
//----------------------------------------------------------------------------
void SynapseGroup::setDefaultMaxConnections()
{
    // If connectivity is procedural
    if(m_MatrixType == SynapseMatrixConnectivity::PROCEDURAL) {
        // If there's no row build code, give an error
        if(m_ConnectivityInitialiser.getSnippet()->getRowBuildCode().empty()) {
            throw std::runtime_error("Cannot use procedural connectivity without specifying connectivity initialisation snippet");
        }

        // If the weight update model has code for postsynaptic-spike triggered updating, give an error
        if(!m_WUModel->getLearnPostCode().empty()) {
            throw std::runtime_error("Procedural connectivity cannot be used for synapse groups with postsynaptic spike-triggered learning");
        }
    }

    // If connectivitity initialisation snippet provides a function to calculate row length, call it
    // **NOTE** only do this for sparse connectivity as this should not be set for bitmasks
    auto calcMaxRowLengthFunc = m_ConnectivityInitialiser.getSnippet()->getCalcMaxRowLengthFunc();
    if(calcMaxRowLengthFunc && (m_MatrixConnectivity == SynapseMatrixConnectivity::SPARSE)) {
        m_MaxConnections = calcMaxRowLengthFunc(getSrcNeuronGroup()->getNumNeurons(), getTrgNeuronGroup()->getNumNeurons(),
                                                m_ConnectivityInitialiser.getParams());
    }
    // Otherwise, default to the size of the target population
    else {
        m_MaxConnections = getTrgNeuronGroup()->getNumNeurons();
    }

    // If connectivitity initialisation snippet provides a function to calculate row length, call it
    // **NOTE** only do this for sparse connectivity as this should not be set for bitmasks
    auto calcMaxColLengthFunc = m_ConnectivityInitialiser.getSnippet()->getCalcMaxColLengthFunc();
    if(calcMaxColLengthFunc && (m_MatrixConnectivity == SynapseMatrixConnectivity::SPARSE)) {
        m_MaxSourceConnections = calcMaxColLengthFunc(getSrcNeuronGroup()->getNumNeurons(), getTrgNeuronGroup()->getNumNeurons(),
                                                      m_ConnectivityInitialiser.getParams());
    }
    // Otherwise, default to the size of the source population
    else {
        m_MaxSourceConnections = getSrcNeuronGroup()->getNumNeurons();
    }

    // Check that the source neuron group supports the desired number of delay steps
    getSrcNeuronGroup()->checkNumDelaySlots(getDelaySteps());
}
