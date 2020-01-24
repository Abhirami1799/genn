#include "code_generator/generateSynapseUpdate.h"

// Standard C++ includes
#include <string>

// GeNN code generator includes
#include "code_generator/codeStream.h"
#include "code_generator/substitutions.h"
#include "code_generator/backendBase.h"
#include "code_generator/groupMerged.h"
#include "code_generator/modelSpecMerged.h"
#include "code_generator/teeStream.h"

//--------------------------------------------------------------------------
// Anonymous namespace
//--------------------------------------------------------------------------
namespace
{
void applySynapseSubstitutions(CodeGenerator::CodeStream &os, std::string code, const std::string &errorContext,
                               const CodeGenerator::SynapseGroupMerged &sg, const CodeGenerator::Substitutions &baseSubs,
                               const ModelSpecInternal &model, const CodeGenerator::BackendBase &backend)
{
    const auto *wu = sg.getArchetype().getWUModel();

    CodeGenerator::Substitutions synapseSubs(&baseSubs);

    // Substitute parameter and derived parameter names
    synapseSubs.addParamValueSubstitution(wu->getParamNames(), sg.getArchetype().getWUParams());
    synapseSubs.addVarValueSubstitution(wu->getDerivedParams(), sg.getArchetype().getWUDerivedParams());
    synapseSubs.addVarNameSubstitution(wu->getExtraGlobalParams(), "", "group.");

    // Substitute names of pre and postsynaptic weight update variables
    const std::string delayedPreIdx = (sg.getArchetype().getDelaySteps() == NO_DELAY) ? synapseSubs["id_pre"] : "preReadDelayOffset + " + baseSubs["id_pre"];
    synapseSubs.addVarNameSubstitution(wu->getPreVars(), "", "group.",
                                       "[" + delayedPreIdx + "]");

    const std::string delayedPostIdx = (sg.getArchetype().getBackPropDelaySteps() == NO_DELAY) ? synapseSubs["id_post"] : "postReadDelayOffset + " + baseSubs["id_post"];
    synapseSubs.addVarNameSubstitution(wu->getPostVars(), "", "group.",
                                       "[" + delayedPostIdx + "]");

    // If weights are individual, substitute variables for values stored in global memory
    if (sg.getArchetype().getMatrixType() & SynapseMatrixWeight::INDIVIDUAL) {
        synapseSubs.addVarNameSubstitution(wu->getVars(), "", "group.",
                                           "[" + synapseSubs["id_syn"] + "]");
    }
    // Otherwise, if weights are procedual
    else if (sg.getArchetype().getMatrixType() & SynapseMatrixWeight::PROCEDURAL) {
        const auto vars = wu->getVars();

        for(size_t k = 0; k < vars.size(); k++) {
            const auto &varInit = sg.getArchetype().getWUVarInitialisers().at(k);

            // If this variable has any initialisation code
            if(!varInit.getSnippet()->getCode().empty()) {
                // Configure variable substitutions
                CodeGenerator::Substitutions varSubs(&synapseSubs);
                varSubs.addVarSubstitution("value", "l" + vars[k].name);
                varSubs.addParamValueSubstitution(varInit.getSnippet()->getParamNames(), varInit.getParams(),
                                                  [k, &sg](size_t p) { return sg.isWUVarInitParamHeterogeneous(k, p); },
                                                  "", "group.", vars[k].name);
                varSubs.addVarValueSubstitution(varInit.getSnippet()->getDerivedParams(), varInit.getDerivedParams(),
                                                [k, &sg](size_t p) { return sg.isWUVarInitDerivedParamHeterogeneous(k, p); },
                                                "", "group.", vars[k].name);

                // Generate variable initialization code
                std::string code = varInit.getSnippet()->getCode();
                varSubs.applyCheckUnreplaced(code, "initVar : merged" + vars[k].name + std::to_string(sg.getIndex()));

                // Declare local variable
                os << vars[k].type << " " << "l" << vars[k].name << ";" << std::endl;

                // Insert code to initialize variable into scope
                {
                    CodeGenerator::CodeStream::Scope b(os);
                    os << code << std::endl;;
                }
            }
        }

        // Substitute variables for newly-declared local variables
        synapseSubs.addVarNameSubstitution(vars, "", "l");
    }
    // Otherwise, substitute variables for constant values
    else {
        synapseSubs.addVarValueSubstitution(wu->getVars(), sg.getArchetype().getWUConstInitVals());
    }

    neuronSubstitutionsInSynapticCode(synapseSubs, sg.getArchetype(), synapseSubs["id_pre"],
                                      synapseSubs["id_post"], model.getDT());

    synapseSubs.apply(code);
    //synapseSubs.applyCheckUnreplaced(code, errorContext + " : " + sg.getName());
    code = CodeGenerator::ensureFtype(code, model.getPrecision());
    os << code;
}
}   // Anonymous namespace

//--------------------------------------------------------------------------
// CodeGenerator
//--------------------------------------------------------------------------
void CodeGenerator::generateSynapseUpdate(CodeStream &os, const MergedEGPMap &mergedEGPs, const ModelSpecMerged &modelMerged, const BackendBase &backend,
                                          bool standaloneModules)
{
    if(standaloneModules) {
        os << "#include \"runner.cc\"" << std::endl;
    }
    else {
        os << "#include \"definitionsInternal.h\"" << std::endl;
    }
    os << "#include \"supportCode.h\"" << std::endl;
    os << std::endl;

    // Generate functions to push merged synapse group structures
    const ModelSpecInternal &model = modelMerged.getModel();
    genMergedGroupPush(os, modelMerged.getMergedSynapseDendriticDelayUpdateGroups(), mergedEGPs, "SynapseDendriticDelayUpdate", backend);
    genMergedGroupPush(os, modelMerged.getMergedPresynapticUpdateGroups(), mergedEGPs, "PresynapticUpdate", backend);
    genMergedGroupPush(os, modelMerged.getMergedPostsynapticUpdateGroups(), mergedEGPs, "PostsynapticUpdate", backend);
    genMergedGroupPush(os, modelMerged.getMergedSynapseDynamicsGroups(), mergedEGPs, "SynapseDynamics", backend);

    // Synaptic update kernels
    backend.genSynapseUpdate(os, modelMerged,
        // Presynaptic weight update threshold
        [&backend, &model](CodeStream &os, const SynapseGroupMerged &sg, Substitutions &baseSubs)
        {
            Substitutions synapseSubs(&baseSubs);

            // Make weight update model substitutions
            synapseSubs.addParamValueSubstitution(sg.getArchetype().getWUModel()->getParamNames(), sg.getArchetype().getWUParams());
            synapseSubs.addVarValueSubstitution(sg.getArchetype().getWUModel()->getDerivedParams(), sg.getArchetype().getWUDerivedParams());
            synapseSubs.addVarNameSubstitution(sg.getArchetype().getWUModel()->getExtraGlobalParams(), "", "group.");

            // Get read offset if required
            const std::string offset = sg.getArchetype().getSrcNeuronGroup()->isDelayRequired() ? "preReadDelayOffset + " : "";
            neuronSubstitutionsInSynapticCode(synapseSubs, sg.getArchetype().getSrcNeuronGroup(), offset, "", baseSubs["id_pre"], "_pre", "Pre");

            // Get event threshold condition code
            std::string code = sg.getArchetype().getWUModel()->getEventThresholdConditionCode();
            synapseSubs.applyCheckUnreplaced(code, "eventThresholdConditionCode");
            code = ensureFtype(code, model.getPrecision());
            os << code;
        },
        // Presynaptic spike
        [&backend, &model](CodeStream &os, const SynapseGroupMerged &sg, Substitutions &baseSubs)
        {
            applySynapseSubstitutions(os, sg.getArchetype().getWUModel()->getSimCode(), "simCode",
                                      sg, baseSubs, model, backend);
        },
        // Presynaptic spike-like event
        [&backend, &model](CodeStream &os, const SynapseGroupMerged &sg, Substitutions &baseSubs)
        {
            applySynapseSubstitutions(os, sg.getArchetype().getWUModel()->getEventCode(), "eventCode",
                                      sg, baseSubs, model, backend);
        },
        // Procedural connectivity
        [&backend, &model](CodeStream &os, const SynapseGroupMerged &sg, Substitutions &baseSubs)
        {
            baseSubs.addFuncSubstitution("endRow", 0, "break");

            // Initialise row building state variables for procedural connectivity
            const auto &connectInit = sg.getArchetype().getConnectivityInitialiser();
            for(const auto &a : connectInit.getSnippet()->getRowBuildStateVars()) {
                os << a.type << " " << a.name << " = " << a.value << ";" << std::endl;
            }

            // Loop through synapses in row
            os << "while(true)";
            {
                CodeStream::Scope b(os);
                Substitutions synSubs(&baseSubs);

                synSubs.addParamValueSubstitution(connectInit.getSnippet()->getParamNames(), connectInit.getParams());
                synSubs.addVarValueSubstitution(connectInit.getSnippet()->getDerivedParams(), connectInit.getDerivedParams());
                synSubs.addVarNameSubstitution(connectInit.getSnippet()->getExtraGlobalParams(), "", "group.");
                synSubs.addVarNameSubstitution(connectInit.getSnippet()->getRowBuildStateVars());

                std::string pCode = connectInit.getSnippet()->getRowBuildCode();
                synSubs.applyCheckUnreplaced(pCode, "proceduralSparseConnectivity : merged " + std::to_string(sg.getIndex()));
                pCode = ensureFtype(pCode, model.getPrecision());

                // Write out code
                os << pCode << std::endl;
            }
        },
        // Postsynaptic learning code
        [&backend, &modelMerged](CodeStream &os, const SynapseGroupMerged &sg, const Substitutions &baseSubs)
        {
            const auto *wum = sg.getArchetype().getWUModel();
            if (!wum->getLearnPostSupportCode().empty()) {
                os << "using namespace " << modelMerged.getPostsynapticUpdateSupportCodeNamespace(wum->getLearnPostSupportCode()) <<  ";" << std::endl;
            }

            applySynapseSubstitutions(os, wum->getLearnPostCode(), "learnPostCode",
                                      sg, baseSubs, modelMerged.getModel(), backend);
        },
        // Synapse dynamics
        [&backend, &modelMerged](CodeStream &os, const SynapseGroupMerged &sg, const Substitutions &baseSubs)
        {
            const auto *wum = sg.getArchetype().getWUModel();
            if (!wum->getSynapseDynamicsSuppportCode().empty()) {
                os << "using namespace " << modelMerged.getSynapseDynamicsSupportCodeNamespace(wum->getSynapseDynamicsSuppportCode()) <<  ";" << std::endl;
            }

            applySynapseSubstitutions(os, wum->getSynapseDynamicsCode(), "synapseDynamics",
                                      sg, baseSubs, modelMerged.getModel(), backend);
        },
        // Push EGP handler
        [&backend, &mergedEGPs](CodeStream &os)
        {
            genScalarEGPPush(os, mergedEGPs, "PresynapticUpdate", backend);
            genScalarEGPPush(os, mergedEGPs, "PostsynapticUpdate", backend);
            genScalarEGPPush(os, mergedEGPs, "SynapseDynamics", backend);
        });
}
