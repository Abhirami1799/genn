//--------------------------------------------------------------------------
/*! \file synapse_support_code_synapse_dynamics/model_new.cc

\brief model definition file that is part of the feature testing
suite of minimal models with known analytic outcomes that are used for continuous integration testing.
*/
//--------------------------------------------------------------------------


#include "modelSpec.h"

//----------------------------------------------------------------------------
// Neuron
//----------------------------------------------------------------------------
class Neuron : public NeuronModels::Base
{
public:
    DECLARE_MODEL(Neuron, 0, 2);

    SET_SIM_CODE("$(x)= $(t)+$(shift);\n");

    SET_VARS({{"x", "scalar"}, {"shift", "scalar"}});
};

IMPLEMENT_MODEL(Neuron);

//----------------------------------------------------------------------------
// WeightUpdateModel
//----------------------------------------------------------------------------
class WeightUpdateModel : public WeightUpdateModels::Base
{
public:
    DECLARE_MODEL(WeightUpdateModel, 0, 1);

    SET_VARS({{"w", "scalar"}});

    SET_SYNAPSE_DYNAMICS_SUPPORT_CODE("__device__ __host__ scalar getWeight(scalar x){ return x; }");
    SET_SYNAPSE_DYNAMICS_CODE("$(w)= getWeight($(x_post));");
};

IMPLEMENT_MODEL(WeightUpdateModel);

void modelDefinition(NNmodel &model)
{
    model.setDT(0.1);
    model.setName("synapse_support_code_synapse_dynamics_new");

    model.addNeuronPopulation<Neuron>("pre", 10, {}, Neuron::VarValues(0.0, 0.0));
    model.addNeuronPopulation<Neuron>("post", 10, {}, Neuron::VarValues(0.0, 0.0));

    string synName= "syn";
    for (int i= 0; i < 10; i++)
    {
        string theName= synName + std::to_string(i);
        model.addSynapsePopulation<WeightUpdateModel, PostsynapticModels::DeltaCurr>(
            theName, SynapseMatrixType::DENSE_INDIVIDUALG, i, "pre", "post",
            {}, WeightUpdateModel::VarValues(0.0),
            {}, {});
    }
    model.setPrecision(GENN_FLOAT);
}
