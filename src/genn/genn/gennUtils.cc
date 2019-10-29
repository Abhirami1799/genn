#include "gennUtils.h"

// Standard C++ includes
#include <algorithm>

// PLOG includes
#include <plog/Log.h>

namespace
{
//--------------------------------------------------------------------------
// GenericFunction
//--------------------------------------------------------------------------
//! Immutable structure for specifying the name and number of
//! arguments of a generic funcion e.g. gennrand_uniform
struct GenericFunction
{
    //! Generic name used to refer to function in user code
    const std::string genericName;

    //! Number of function arguments
    const unsigned int numArguments;
};


GenericFunction randomFuncs[] = {
    {"gennrand_uniform", 0},
    {"gennrand_normal", 0},
    {"gennrand_exponential", 0},
    {"gennrand_log_normal", 2},
    {"gennrand_gamma", 1}
};
}

//--------------------------------------------------------------------------
// Utils
//--------------------------------------------------------------------------
namespace Utils
{
bool isRNGRequired(const std::string &code)
{
    // Loop through random functions
    for(const auto &r : randomFuncs) {
        // If this function takes no arguments, return true if
        // generic function name enclosed in $() markers is found
        if(r.numArguments == 0) {
            if(code.find("$(" + r.genericName + ")") != std::string::npos) {
                return true;
            }
        }
        // Otherwise, return true if generic function name
        // prefixed by $( and suffixed with comma is found
        else if(code.find("$(" + r.genericName + ",") != std::string::npos) {
            return true;
        }
    }
    return false;

}
//--------------------------------------------------------------------------
bool isInitRNGRequired(const std::vector<Models::VarInit> &varInitialisers)
{
    // Return true if any of these variable initialisers require an RNG
    return std::any_of(varInitialisers.cbegin(), varInitialisers.cend(),
                       [](const Models::VarInit &varInit)
                       {
                           return isRNGRequired(varInit.getSnippet()->getCode());
                       });
}
//--------------------------------------------------------------------------
void initialiseLegacyImplementation(const std::vector<double> &params, const std::vector<Models::VarInit> &initialisers,
                                    std::vector<Models::VarInit> &combinedInitialisers, std::vector<VarImplementation> &implementation,
                                    VarImplementation defaultVarImplementation)
{
    // Reserve combined initialisers and implementations to be large enough for COMBINED vars
    combinedInitialisers.reserve(params.size() + initialisers.size());
    implementation.reserve(params.size() + initialisers.size());

    // Transform parameter values into constant var initialisers
    std::transform(params.cbegin(), params.cend(), std::back_inserter(combinedInitialisers),
                   [](double v){ return Models::VarInit(v); });

    // Copy variable initialisers after
    combinedInitialisers.insert(combinedInitialisers.end(), initialisers.cbegin(), initialisers.cend());

    // Implement all parameters as GLOBAL
    implementation.insert(implementation.end(), params.size(), VarImplementation::GLOBAL);

    // Implement all variables using default implementation
    implementation.insert(implementation.end(), initialisers.size(), defaultVarImplementation);
}
//--------------------------------------------------------------------------
void autoDetermineImplementation(const std::vector<Models::VarInit> &initialisers, const Models::Base::VarVec &vars,
                                 std::vector<VarImplementation> &implementations)
{
    // Assert that all sizes match
    assert(initialisers.size() == vars.size());

    // Reserve implementations arra
    implementations.reserve(initialisers.size());

    // Implement varaibles that are read only and constant as GLOBAL and others as individual
    std::transform(initialisers.cbegin(), initialisers.cend(), vars.cbegin(), std::back_inserter(implementations),
                   [](const Models::VarInit &varInit, const Models::Base::Var &var)
                   {
                       if(var.access == VarAccess::READ_ONLY && varInit.isConstant()) {
                           LOGD << "Defaulting variable '" << var.name << "' to GLOBAL implementation";
                           return VarImplementation::GLOBAL;
                       }
                       else {
                           LOGD << "Defaulting variable '" << var.name << "' to INDIVIDUAL implementation";
                           return VarImplementation::INDIVIDUAL;
                       }
                   });
}
//--------------------------------------------------------------------------
void calcDerivedParamVal(const Models::Base *model, const std::vector<Models::VarInit> &initialisers,
                         const std::vector<VarImplementation> &implementation, double dt,
                         std::vector<double> &derivedParamValues)
{
    auto paramNames = model->getParamNames();
    auto derivedParams = model->getDerivedParams();
    auto derivedParamsNamed = model->getDerivedParamsNamed();
    auto combinedVars =  model->getCombinedVars();

    // Reserve vector to hold derived parameters
    derivedParamValues.reserve(derivedParams.size() + derivedParamsNamed.size());

    // If there are any legacy parameter names
    if(!paramNames.empty()) {
        // initialiseLegacyImplementation should have added them to beginning of var initialisers
        std::vector<double> parValues;
        parValues.reserve(paramNames.size());
        std::transform(paramNames.cbegin(), paramNames.cend(), initialisers.cbegin(), std::back_inserter(parValues),
                       [](const std::string&, const Models::VarInit &v){ return v.getConstantValue(); });

        // Loop through derived parameters
        for(const auto &d : derivedParams) {
            derivedParamValues.push_back(d.func(parValues, dt));
        }
    }

    // Loop through variables and their initializers and build dictionary of the
    // constant value associated with variables implemented as GLOBAL
    std::map<std::string, double> varValues;
    auto var = combinedVars.cbegin();
    auto varInit = initialisers.cbegin();
    auto varImpl = implementation.cbegin();
    for(; var != combinedVars.cend(); var++, varInit++, varImpl++) {
        if(*varImpl == VarImplementation::GLOBAL) {
            varValues.emplace(var->name, varInit->getConstantValue());
        }
    }

    // Loop through named derived parameters
    for(const auto &d : derivedParamsNamed) {
        derivedParamValues.push_back(d.func(varValues, dt));
    }

}
//--------------------------------------------------------------------------
bool isTypePointer(const std::string &type)
{
    return (type.back() == '*');
}
//--------------------------------------------------------------------------
std::string getUnderlyingType(const std::string &type)
{
    // Check that type is a pointer type
    assert(isTypePointer(type));

    // Return string without last character
    return type.substr(0, type.length() - 1);
}
}   // namespace utils
