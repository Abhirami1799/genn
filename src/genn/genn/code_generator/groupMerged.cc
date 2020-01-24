#include "code_generator/groupMerged.h"

// PLOG includes
#include <plog/Log.h>

// GeNN includes
#include "modelSpecInternal.h"

//----------------------------------------------------------------------------
// CodeGenerator::NeuronGroupMerged
//----------------------------------------------------------------------------
CodeGenerator::NeuronGroupMerged::NeuronGroupMerged(size_t index, bool init, const std::vector<std::reference_wrapper<const NeuronGroupInternal>> &groups)
:   CodeGenerator::GroupMerged<NeuronGroupInternal>(index, groups)
{
    // Build vector of vectors containing each child group's merged in syns, ordered to match those of the archetype group
    orderNeuronGroupChildren(m_SortedMergedInSyns, &NeuronGroupInternal::getMergedInSyn,
                             [init](const std::pair<SynapseGroupInternal *, std::vector<SynapseGroupInternal *>> &a,
                                    const std::pair<SynapseGroupInternal *, std::vector<SynapseGroupInternal *>> &b)
                             {
                                 return init ? a.first->canPSInitBeMerged(*b.first) : a.first->canPSBeMerged(*b.first);
                             });

    // Build vector of vectors containing each child group's current sources, ordered to match those of the archetype group
    orderNeuronGroupChildren(m_SortedCurrentSources, &NeuronGroupInternal::getCurrentSources,
                             [init](const CurrentSourceInternal *a, const CurrentSourceInternal *b)
                             {
                                 return init ? a->canInitBeMerged(*b) : a->canBeMerged(*b);
                             });

    // Build vector of vectors containing each child group's incoming synapse groups
    // with postsynaptic updates, ordered to match those of the archetype group
    const auto inSynWithPostCode = getArchetype().getInSynWithPostCode();
    orderNeuronGroupChildren(inSynWithPostCode, m_SortedInSynWithPostCode, &NeuronGroupInternal::getInSynWithPostCode,
                             [init](const SynapseGroupInternal *a, const SynapseGroupInternal *b)
                             {
                                 return init ? a->canWUPostInitBeMerged(*b) : a->canWUPostBeMerged(*b);
                             });

    // Build vector of vectors containing each child group's incoming synapse groups
    // with postsynaptic updates, ordered to match those of the archetype group
    const auto outSynWithPreCode = getArchetype().getOutSynWithPreCode();
    orderNeuronGroupChildren(outSynWithPreCode, m_SortedOutSynWithPreCode, &NeuronGroupInternal::getOutSynWithPreCode,
                             [init](const SynapseGroupInternal *a, const SynapseGroupInternal *b)
                             {
                                 return init ? a->canWUPreInitBeMerged(*b) : a->canWUPreBeMerged(*b);
                             });
}
//----------------------------------------------------------------------------
std::string CodeGenerator::NeuronGroupMerged::getCurrentQueueOffset() const
{
    assert(getArchetype().isDelayRequired());
    return "(*group.spkQuePtr * group.numNeurons)";
}
//----------------------------------------------------------------------------
std::string CodeGenerator::NeuronGroupMerged::getPrevQueueOffset() const
{
    assert(getArchetype().isDelayRequired());
    return "(((*group.spkQuePtr + " + std::to_string(getArchetype().getNumDelaySlots() - 1) + ") % " + std::to_string(getArchetype().getNumDelaySlots()) + ") * group.numNeurons)";
}
//----------------------------------------------------------------------------
bool CodeGenerator::NeuronGroupMerged::isParamHeterogeneous(size_t index) const
{
    return CodeGenerator::GroupMerged<NeuronGroupInternal>::isParamValueHeterogeneous(
        index, [](const NeuronGroupInternal &ng) { return ng.getParams(); });
}
//----------------------------------------------------------------------------
bool CodeGenerator::NeuronGroupMerged::isDerivedParamHeterogeneous(size_t index) const
{
    return CodeGenerator::GroupMerged<NeuronGroupInternal>::isParamValueHeterogeneous(
        index, [](const NeuronGroupInternal &ng) { return ng.getDerivedParams(); });
}
//----------------------------------------------------------------------------
bool CodeGenerator::NeuronGroupMerged::isCurrentSourceParamHeterogeneous(size_t childIndex, size_t paramIndex) const
{
    // If parameter isn't referenced in code, there's no point implementing it hetereogeneously!
    const auto *csm = m_SortedCurrentSources[0][childIndex]->getCurrentSourceModel();
    const std::string paramName = csm->getParamNames().at(paramIndex);
    if(csm->getInjectionCode().find("$(" + paramName + ")") == std::string::npos) {
        return false;
    }
    // Otherwise, return whether values across all groups are heterogeneous
    else {
        return isChildParamValueHeterogeneous(childIndex, paramIndex, m_SortedCurrentSources,
                                              &CurrentSourceInternal::getParams);
    }
}
//----------------------------------------------------------------------------
bool CodeGenerator::NeuronGroupMerged::isCurrentSourceDerivedParamHeterogeneous(size_t childIndex, size_t paramIndex) const
{
    // If derived parameter isn't referenced in code, there's no point implementing it hetereogeneously!
    const auto *csm = m_SortedCurrentSources[0][childIndex]->getCurrentSourceModel();
    const std::string derivedParamName = csm->getDerivedParams().at(paramIndex).name;
    if(csm->getInjectionCode().find("$(" + derivedParamName + ")") == std::string::npos) {
        return false;
    }
    // Otherwise, return whether values across all groups are heterogeneous
    else {
        return isChildParamValueHeterogeneous(childIndex, paramIndex, m_SortedCurrentSources, 
                                              &CurrentSourceInternal::getDerivedParams);
    }
}
//----------------------------------------------------------------------------
bool CodeGenerator::NeuronGroupMerged::isInSynWithPostCodeParamHeterogeneous(size_t childIndex, size_t paramIndex) const
{
    // If parameter isn't referenced in code, there's no point implementing it hetereogeneously!
    const auto *wum = m_SortedInSynWithPostCode[0][childIndex]->getWUModel();
    const std::string paramName = wum->getParamNames().at(paramIndex);
    if(wum->getPostSpikeCode().find("$(" + paramName + ")") == std::string::npos) {
        return false;
    }
    // Otherwise, return whether values across all groups are heterogeneous
    else {
        return isChildParamValueHeterogeneous(childIndex, paramIndex, m_SortedInSynWithPostCode,
                                              &SynapseGroupInternal::getWUParams);
    }
}
//----------------------------------------------------------------------------
bool CodeGenerator::NeuronGroupMerged::isInSynWithPostCodeDerivedParamHeterogeneous(size_t childIndex, size_t paramIndex) const
{
    // If parameter isn't referenced in code, there's no point implementing it hetereogeneously!
    const auto *wum = m_SortedInSynWithPostCode[0][childIndex]->getWUModel();
    const std::string derivedParamName = wum->getDerivedParams().at(paramIndex).name;
    if(wum->getPostSpikeCode().find("$(" + derivedParamName + ")") == std::string::npos) {
        return false;
    }
    // Otherwise, return whether values across all groups are heterogeneous
    else {
        return isChildParamValueHeterogeneous(childIndex, paramIndex, m_SortedInSynWithPostCode,
                                              &SynapseGroupInternal::getWUDerivedParams);
    }
}
//----------------------------------------------------------------------------
bool CodeGenerator::NeuronGroupMerged::isOutSynWithPreCodeParamHeterogeneous(size_t childIndex, size_t paramIndex) const
{
    // If parameter isn't referenced in code, there's no point implementing it hetereogeneously!
    const auto *wum = m_SortedOutSynWithPreCode[0][childIndex]->getWUModel();
    const std::string paramName = wum->getParamNames().at(paramIndex);
    if(wum->getPreSpikeCode().find("$(" + paramName + ")") == std::string::npos) {
        return false;
    }
    // Otherwise, return whether values across all groups are heterogeneous
    else {
        return isChildParamValueHeterogeneous(childIndex, paramIndex, m_SortedOutSynWithPreCode,
                                              &SynapseGroupInternal::getWUParams);
    }
}
//----------------------------------------------------------------------------
bool CodeGenerator::NeuronGroupMerged::isOutSynWithPreCodeDerivedParamHeterogeneous(size_t childIndex, size_t paramIndex) const
{
    // If parameter isn't referenced in code, there's no point implementing it hetereogeneously!
    const auto *wum = m_SortedOutSynWithPreCode[0][childIndex]->getWUModel();
    const std::string derivedParamName = wum->getDerivedParams().at(paramIndex).name;
    if(wum->getPreSpikeCode().find("$(" + derivedParamName + ")") == std::string::npos) {
        return false;
    }
    // Otherwise, return whether values across all groups are heterogeneous
    else {
        return isChildParamValueHeterogeneous(childIndex, paramIndex, m_SortedOutSynWithPreCode,
                                              &SynapseGroupInternal::getWUDerivedParams);
    }
}

//----------------------------------------------------------------------------
// CodeGenerator::SynapseGroupMerged
//----------------------------------------------------------------------------
std::string CodeGenerator::SynapseGroupMerged::getPresynapticAxonalDelaySlot() const
{
    assert(getArchetype().getSrcNeuronGroup()->isDelayRequired());

    const unsigned int numDelaySteps = getArchetype().getDelaySteps();
    if(numDelaySteps == 0) {
        return "(*group.srcSpkQuePtr)";
    }
    else {
        const unsigned int numSrcDelaySlots = getArchetype().getSrcNeuronGroup()->getNumDelaySlots();
        return "((*group.srcSpkQuePtr + " + std::to_string(numSrcDelaySlots - numDelaySteps) + ") % " + std::to_string(numSrcDelaySlots) + ")";
    }
}
//----------------------------------------------------------------------------
std::string CodeGenerator::SynapseGroupMerged::getPostsynapticBackPropDelaySlot() const
{
    assert(getArchetype().getTrgNeuronGroup()->isDelayRequired());

    const unsigned int numBackPropDelaySteps = getArchetype().getBackPropDelaySteps();
    if(numBackPropDelaySteps == 0) {
        return "(*group.trgSpkQuePtr)";
    }
    else {
        const unsigned int numTrgDelaySlots = getArchetype().getTrgNeuronGroup()->getNumDelaySlots();
        return "((*group.trgSpkQuePtr + " + std::to_string(numTrgDelaySlots - numBackPropDelaySteps) + ") % " + std::to_string(numTrgDelaySlots) + ")";
    }
}
//----------------------------------------------------------------------------
std::string CodeGenerator::SynapseGroupMerged::getDendriticDelayOffset(const std::string &offset) const
{
    assert(getArchetype().isDendriticDelayRequired());

    if(offset.empty()) {
        return "(*group.denDelayPtr * group.numTrgNeurons) + ";
    }
    else {
        return "(((*group.denDelayPtr + " + offset + ") % " + std::to_string(getArchetype().getMaxDendriticDelayTimesteps()) + ") * group.numTrgNeurons) + ";
    }
}