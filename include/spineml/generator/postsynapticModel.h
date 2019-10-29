#pragma once

#pragma once

// Standard includes
#include <set>
#include <string>

// GeNN includes
#include "postsynapticModels.h"

// Spine ML generator includes
#include "modelCommon.h"

// Forward declarations
namespace SpineMLGenerator
{
    class NeuronModel;
    class WeightUpdateModel;

    namespace ModelParams
    {
        class Postsynaptic;
    }
}

namespace pugi
{
    class xml_node;
}

//----------------------------------------------------------------------------
// SpineMLGenerator::PostsynapticModel
//----------------------------------------------------------------------------
namespace SpineMLGenerator
{
class PostsynapticModel : public PostsynapticModels::Base
{
public:
    PostsynapticModel(const ModelParams::Postsynaptic &params,
                      const pugi::xml_node &componentClass,
                      const NeuronModel *trgNeuronModel,
                      const WeightUpdateModel *weightUpdateModel);

    //------------------------------------------------------------------------
    // Typedefines
    //------------------------------------------------------------------------
    typedef SpineMLGenerator::VarValues<PostsynapticModel> VarValues;

    //------------------------------------------------------------------------
    // PostsynapticModels::Base virtuals
    //------------------------------------------------------------------------
    virtual std::string getDecayCode() const override{ return m_DecayCode; }
    virtual std::string getApplyInputCode() const override{ return m_ApplyInputCode; }
    virtual Models::Base::VarVec getVars() const override{ return m_Vars; }
    virtual Models::Base::DerivedParamNamedVec getDerivedParamsNamed() const override{ return m_DerivedParams; }

    //------------------------------------------------------------------------
    // Public API
    //------------------------------------------------------------------------
    unsigned int getInitialRegimeID() const
    {
        return m_InitialRegimeID;
    }

    //------------------------------------------------------------------------
    // Constants
    //------------------------------------------------------------------------
    static const char *componentClassName;

private:
    //------------------------------------------------------------------------
    // Members
    //------------------------------------------------------------------------
    std::string m_DecayCode;
    std::string m_ApplyInputCode;

    Models::Base::VarVec m_Vars;
    Models::Base::DerivedParamNamedVec m_DerivedParams;

    unsigned int m_InitialRegimeID;
};
}   // namespace SpineMLGenerator
