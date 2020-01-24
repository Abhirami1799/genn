#pragma once

// Standard includes
#include <functional>
#include <vector>

// GeNN includes
#include "gennExport.h"
#include "currentSourceInternal.h"
#include "neuronGroupInternal.h"
#include "synapseGroupInternal.h"

//----------------------------------------------------------------------------
// CodeGenerator::GroupMerged
//----------------------------------------------------------------------------
//! Very thin wrapper around a number of groups which have been merged together
namespace CodeGenerator
{
template<typename G>
class GroupMerged
{
public:
    //------------------------------------------------------------------------
    // Typedefines
    //------------------------------------------------------------------------
    typedef G GroupInternal;

    GroupMerged(size_t index, const std::vector<std::reference_wrapper<const GroupInternal>> &groups)
    :   m_Index(index), m_Groups(groups)
    {}

    //------------------------------------------------------------------------
    // Public API
    //------------------------------------------------------------------------
    size_t getIndex() const { return m_Index; }

    //! Get 'archetype' neuron group - it's properties represent those of all other merged neuron groups
    const GroupInternal &getArchetype() const { return m_Groups.front().get(); }

    //! Gets access to underlying vector of neuron groups which have been merged
    const std::vector<std::reference_wrapper<const GroupInternal>> &getGroups() const{ return m_Groups; }

protected:
    //------------------------------------------------------------------------
    // Protected methods
    //------------------------------------------------------------------------
    //! Helper to test whether parameter values are heterogeneous within merged group
    template<typename P>
    bool isParamValueHeterogeneous(size_t index, P getParamValuesFn) const
    {
        // Get value of parameter in archetype group
        const double archetypeValue = getParamValuesFn(getArchetype()).at(index);

        // Return true if any parameter values differ from the archetype value
        return std::any_of(getGroups().cbegin(), getGroups().cend(),
                           [archetypeValue, index, getParamValuesFn](const GroupInternal &g)
                           {
                               return (getParamValuesFn(g).at(index) != archetypeValue);
                           });
    }

    //! Helper to test whether a variable initialiser param values are heterogeneous across merged group
    template<typename V>
    bool isVarInitParamHeterogeneous(size_t varIndex, size_t paramIndex, 
                                     V getVarInitialisers) const
    {
        // If parameter isn't referenced in code, there's no point implementing it hetereogeneously!
        const auto *varInitSnippet = (getArchetype().*getVarInitialisers)().at(varIndex).getSnippet();
        const std::string paramName = varInitSnippet->getParamNames().at(paramIndex);
        if(varInitSnippet->getCode().find("$(" + paramName + ")") == std::string::npos) {
            return false;
        }
        // Otherwise, return whether values across all groups are heterogeneous
        else {
            return isParamValueHeterogeneous(
                paramIndex,
                [varIndex, getVarInitialisers](const G &sg)
                {
                    return (sg.*getVarInitialisers)().at(varIndex).getParams();
                });
        }
    }

    //! Helper to test whether a variable initialiser derived param values are heterogeneous across merged group
    template<typename V>
    bool isVarInitDerivedParamHeterogeneous(size_t varIndex, size_t paramIndex,
                                            V getVarInitialisers) const
    {
        // If parameter isn't referenced in code, there's no point implementing it hetereogeneously!
        const auto *varInitSnippet = (getArchetype().*getVarInitialisers)().at(varIndex).getSnippet();
        const std::string derivedParamName = varInitSnippet->getDerivedParams().at(paramIndex).name;
        if(varInitSnippet->getCode().find("$(" + derivedParamName + ")") == std::string::npos) {
            return false;
        }
        // Otherwise, return whether values across all groups are heterogeneous
        else {
            return isParamValueHeterogeneous(
                paramIndex,
                [varIndex, getVarInitialisers](const G &sg)
                {
                    return (sg.*getVarInitialisers)().at(varIndex).getDerivedParams();
                });
        }
    }

private:
    //------------------------------------------------------------------------
    // Members
    //------------------------------------------------------------------------
    const size_t m_Index;
    std::vector<std::reference_wrapper<const GroupInternal>> m_Groups;
};

//----------------------------------------------------------------------------
// CodeGenerator::NeuronGroupMerged
//----------------------------------------------------------------------------
class GENN_EXPORT NeuronGroupMerged : public GroupMerged<NeuronGroupInternal>
{
public:
    NeuronGroupMerged(size_t index, bool init, const std::vector<std::reference_wrapper<const NeuronGroupInternal>> &groups);

    //------------------------------------------------------------------------
    // Public API
    //------------------------------------------------------------------------
    //! Get the expression to calculate the queue offset for accessing state of variables this timestep
    std::string getCurrentQueueOffset() const;

    //! Get the expression to calculate the queue offset for accessing state of variables in previous timestep
    std::string getPrevQueueOffset() const;

    //! Is the parameter implemented as a heterogeneous parameter?
    bool isParamHeterogeneous(size_t index) const;

    //! Is the derived parameter implemented as a heterogeneous parameter?
    bool isDerivedParamHeterogeneous(size_t index) const;

    //! Is the neuron variable initialization parameter implemented as a heterogeneous parameter?
    bool isVarInitParamHeterogeneous(size_t varIndex, size_t paramIndex) const
    {
        return GroupMerged<NeuronGroupInternal>::isVarInitParamHeterogeneous(
            varIndex, paramIndex, &NeuronGroupInternal::getVarInitialisers);
    }

    //! Is the neuron variable initialization derived parameter implemented as a heterogeneous parameter?
    bool isVarInitDerivedParamHeterogeneous(size_t varIndex, size_t paramIndex) const
    {
        return GroupMerged<NeuronGroupInternal>::isVarInitDerivedParamHeterogeneous(
            varIndex, paramIndex, &NeuronGroupInternal::getVarInitialisers);
    }

    //! Is the current source parameter implemented as a heterogeneous parameter?
    bool isCurrentSourceParamHeterogeneous(size_t childIndex, size_t paramIndex) const;

    //! Is the current source derived parameter implemented as a heterogeneous parameter?
    bool isCurrentSourceDerivedParamHeterogeneous(size_t childIndex, size_t paramIndex) const;

    //! Is the current source variable initialisation parameter implemented as a heterogeneous parameter?
    bool isCurrentSourceVarInitParamHeterogeneous(size_t childIndex, size_t varIndex, size_t paramIndex) const
    {
        return isChildVarInitParamHeterogeneous(childIndex, varIndex, paramIndex, 
                                                m_SortedCurrentSources, &CurrentSourceInternal::getVarInitialisers);
    }

    //! Is the current source variable initialisation derived parameter implemented as a heterogeneous parameter?
    bool isCurrentSourceVarInitDerivedParamHeterogeneous(size_t childIndex, size_t varIndex, size_t paramIndex) const
    {
        return isChildVarInitDerivedParamHeterogeneous(childIndex, varIndex, paramIndex,
                                                       m_SortedCurrentSources, &CurrentSourceInternal::getVarInitialisers);
    }

    //! Is the insyn postsynaptic parameter implemented as a heterogeneous parameter?
    bool isInSynWithPostCodeParamHeterogeneous(size_t childIndex, size_t paramIndex) const;

    //! Is the insyn postsynaptic derived parameter implemented as a heterogeneous parameter?
    bool isInSynWithPostCodeDerivedParamHeterogeneous(size_t childIndex, size_t paramIndex) const;

    //! Is the insyn postsynaptic variable initialisation parameter implemented as a heterogeneous parameter?
    bool isInSynWithPostCodeVarInitParamHeterogeneous(size_t childIndex, size_t varIndex, size_t paramIndex) const
    {
        return isChildVarInitParamHeterogeneous(childIndex, varIndex, paramIndex,
                                                m_SortedInSynWithPostCode, &SynapseGroupInternal::getWUPostVarInitialisers);
    }

    //! Is the insyn postsynaptic variable initialisation derived parameter implemented as a heterogeneous parameter?
    bool isInSynWithPostCodeVarInitDerivedParamHeterogeneous(size_t childIndex, size_t varIndex, size_t paramIndex) const
    {
        return isChildVarInitDerivedParamHeterogeneous(childIndex, varIndex, paramIndex,
                                                       m_SortedInSynWithPostCode, &SynapseGroupInternal::getWUPostVarInitialisers);
    }

    //! Is the outsyn presynaptic parameter implemented as a heterogeneous parameter?
    bool isOutSynWithPreCodeParamHeterogeneous(size_t childIndex, size_t paramIndex) const;

    //! Is the outsyn presynaptic derived parameter implemented as a heterogeneous parameter?
    bool isOutSynWithPreCodeDerivedParamHeterogeneous(size_t childIndex, size_t paramIndex) const;
   
    //! Is the outsyn presynaptic variable initialisation parameter implemented as a heterogeneous parameter?
    bool isOutSynWithPreCodeVarInitParamHeterogeneous(size_t childIndex, size_t varIndex, size_t paramIndex) const
    {
        return isChildVarInitParamHeterogeneous(childIndex, varIndex, paramIndex,
                                                m_SortedOutSynWithPreCode, &SynapseGroupInternal::getWUPreVarInitialisers);
    }

    //! Is the outsyn presynaptic variable initialisation derived parameter implemented as a heterogeneous parameter?
    bool isOutSynWithPreCodeVarInitDerivedParamHeterogeneous(size_t childIndex, size_t varIndex, size_t paramIndex) const
    {
        return isChildVarInitDerivedParamHeterogeneous(childIndex, varIndex, paramIndex,
                                                       m_SortedOutSynWithPreCode, &SynapseGroupInternal::getWUPreVarInitialisers);
    }

    const std::vector<std::vector<std::pair<SynapseGroupInternal *, std::vector<SynapseGroupInternal *>>>> &getSortedMergedInSyns() const{ return m_SortedMergedInSyns; }
    const std::vector<std::vector<CurrentSourceInternal *>> &getSortedCurrentSources() const { return m_SortedCurrentSources; }
    const std::vector<std::vector<SynapseGroupInternal *>> &getSortedInSynWithPostCode() const { return m_SortedInSynWithPostCode; }
    const std::vector<std::vector<SynapseGroupInternal *>> &getSortedOutSynWithPreCode() const{ return m_SortedOutSynWithPreCode; }

private:
    //------------------------------------------------------------------------
    // Private methods
    //------------------------------------------------------------------------
    template<typename T, typename G, typename C>
    void orderNeuronGroupChildren(const std::vector<T> &archetypeChildren,
                                  std::vector<std::vector<T>> &sortedGroupChildren,
                                  G getVectorFunc, C isCompatibleFunc) const
    {
        // Reserve vector of vectors to hold children for all neuron groups, in archetype order
        sortedGroupChildren.reserve(archetypeChildren.size());

        // Loop through groups
        for(const auto &g : getGroups()) {
            // Make temporary copy of this group's children
            std::vector<T> tempChildren((g.get().*getVectorFunc)());

            assert(tempChildren.size() == archetypeChildren.size());

            // Reserve vector for this group's children
            sortedGroupChildren.emplace_back();
            sortedGroupChildren.back().reserve(tempChildren.size());

            // Loop through archetype group's children
            for(const auto &archetypeG : archetypeChildren) {
                // Find compatible child in temporary list
                const auto otherChild = std::find_if(tempChildren.cbegin(), tempChildren.cend(),
                                                     [archetypeG, isCompatibleFunc](const T &g)
                                                     {
                                                         return isCompatibleFunc(archetypeG, g);
                                                     });
                assert(otherChild != tempChildren.cend());

                // Add pointer to vector of compatible merged in syns
                sortedGroupChildren.back().push_back(*otherChild);

                // Remove from original vector
                tempChildren.erase(otherChild);
            }
        }
    }
    
    template<typename T, typename G, typename C>
    void orderNeuronGroupChildren(std::vector<std::vector<T>> &sortedGroupChildren,
                                  G getVectorFunc, C isCompatibleFunc) const
    {
        const std::vector<T> &archetypeChildren = (getArchetype().*getVectorFunc)();
        orderNeuronGroupChildren(archetypeChildren, sortedGroupChildren, getVectorFunc, isCompatibleFunc);
    }
    
    template<typename T, typename G>
    bool isChildParamValueHeterogeneous(size_t childIndex, size_t paramIndex, const std::vector<std::vector<T>> &sortedGroupChildren,
                                        G getParamValuesFn) const
    {
        // Get value of archetype derived parameter
        const double firstValue = (sortedGroupChildren[0][childIndex]->*getParamValuesFn)().at(paramIndex);

        // Loop through groups within merged group
        for(size_t i = 0; i < sortedGroupChildren.size(); i++) {
            const auto group = sortedGroupChildren[i][childIndex];
            if((group->*getParamValuesFn)().at(paramIndex) != firstValue) {
                return true;
            }
        }
        return false;
    }
    
    template<typename T, typename G>
    bool isChildVarInitParamHeterogeneous(size_t childIndex, size_t varIndex, size_t paramIndex, 
                                          const std::vector<std::vector<T>> &sortedGroupChildren,
                                          G getVarInitialisersFn) const
    {
        // If parameter isn't referenced in code, there's no point implementing it hetereogeneously!
        const auto &firstVarInit = (sortedGroupChildren[0][childIndex]->*getVarInitialisersFn)().at(varIndex);
        const std::string paramName = firstVarInit.getSnippet()->getParamNames().at(paramIndex);
        if(firstVarInit.getSnippet()->getCode().find("$(" + paramName + ")") == std::string::npos) {
            return false;
        }
        // Otherwise, return whether values across all groups are heterogeneous
        else {
            // Loop through groups within merged group
            const double firstValue = firstVarInit.getParams().at(paramIndex);
            for(size_t i = 0; i < sortedGroupChildren.size(); i++) {
                const auto &varInit = (sortedGroupChildren[i][childIndex]->*getVarInitialisersFn)().at(varIndex);
                if(varInit.getParams().at(paramIndex) != firstValue) {
                    return true;
                }
            }
            return false;
        }
    }

    template<typename T, typename G>
    bool isChildVarInitDerivedParamHeterogeneous(size_t childIndex, size_t varIndex, size_t paramIndex,
                                                 const std::vector<std::vector<T>> &sortedGroupChildren,
                                                 G getVarInitialisersFn) const
    {
        // If parameter isn't referenced in code, there's no point implementing it hetereogeneously!
        const auto &firstVarInit = (sortedGroupChildren[0][childIndex]->*getVarInitialisersFn)().at(varIndex);
        const std::string derivedParamName = firstVarInit.getSnippet()->getDerivedParams().at(paramIndex).name;
        if(firstVarInit.getSnippet()->getCode().find("$(" + derivedParamName + ")") == std::string::npos) {
            return false;
        }
        // Otherwise, return whether values across all groups are heterogeneous
        else {
            // Loop through groups within merged group
            const double firstValue = firstVarInit.getDerivedParams().at(paramIndex);
            for(size_t i = 0; i < sortedGroupChildren.size(); i++) {
                const auto &varInit = (sortedGroupChildren[i][childIndex]->*getVarInitialisersFn)().at(varIndex);
                if(varInit.getDerivedParams().at(paramIndex) != firstValue) {
                    return true;
                }
            }
            return false;
        }
    }

    //------------------------------------------------------------------------
    // Members
    //------------------------------------------------------------------------
    std::vector<std::vector<std::pair<SynapseGroupInternal *, std::vector<SynapseGroupInternal *>>>> m_SortedMergedInSyns;
    std::vector<std::vector<CurrentSourceInternal*>> m_SortedCurrentSources;
    std::vector<std::vector<SynapseGroupInternal *>> m_SortedInSynWithPostCode;
    std::vector<std::vector<SynapseGroupInternal *>> m_SortedOutSynWithPreCode;
};

//----------------------------------------------------------------------------
// CodeGenerator::SynapseGroupMerged
//----------------------------------------------------------------------------
class GENN_EXPORT SynapseGroupMerged : public GroupMerged<SynapseGroupInternal>
{
public:
    SynapseGroupMerged(size_t index, bool, const std::vector<std::reference_wrapper<const SynapseGroupInternal>> &groups)
    :   GroupMerged<SynapseGroupInternal>(index, groups)
    {}

    //------------------------------------------------------------------------
    // Public API
    //------------------------------------------------------------------------
    //! Get the expression to calculate the delay slot for accessing
    //! Presynaptic neuron state variables, taking into account axonal delay
    std::string getPresynapticAxonalDelaySlot() const;

    //! Get the expression to calculate the delay slot for accessing
    //! Postsynaptic neuron state variables, taking into account back propagation delay
    std::string getPostsynapticBackPropDelaySlot() const;

    std::string getDendriticDelayOffset(const std::string &offset = "") const;

    //! Is the weight update model variable initialization parameter implemented as a heterogeneous parameter?
    bool isWUVarInitParamHeterogeneous(size_t varIndex, size_t paramIndex) const
    {
        return isVarInitParamHeterogeneous(varIndex, paramIndex, &SynapseGroupInternal::getWUVarInitialisers);
    }
    
    //! Is the weight update model variable initialization derived parameter implemented as a heterogeneous parameter?
    bool isWUVarInitDerivedParamHeterogeneous(size_t varIndex, size_t paramIndex) const
    {
        return isVarInitDerivedParamHeterogeneous(varIndex, paramIndex, &SynapseGroupInternal::getWUVarInitialisers);
    }
};
}   // namespace CodeGenerator
