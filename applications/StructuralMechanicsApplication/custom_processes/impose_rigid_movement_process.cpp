// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "utilities/openmp_utils.h"
#include "custom_processes/impose_rigid_movement_process.h"

namespace Kratos
{
ImposeRigidMovementProcess::ImposeRigidMovementProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"             : "please_specify_model_part_name",
        "new_model_part_name"         : "Rigid_Movement_ModelPart",
        "master_variable_name"        : "DISPLACEMENT",
        "slave_variable_name"         : "",
        "relation"                    : 1.0,
        "constant"                    : 0.0,
        "master_node_id"              : 0
    })" );

    mThisParameters.ValidateAndAssignDefaults(default_parameters);
}

/***********************************************************************************/
/***********************************************************************************/

void ImposeRigidMovementProcess::Execute()
{
    KRATOS_TRY

    // We execute the different steps
    ExecuteInitialize();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void ImposeRigidMovementProcess::ExecuteInitialize()
{
    KRATOS_TRY

    // Getting model parts
    ModelPart& r_root_model_part = mrThisModelPart.GetRootModelPart();
    Model& r_model = r_root_model_part.GetModel();
    ModelPart& r_model_part = r_model.GetModelPart(mThisParameters["model_part_name"].GetString());
    const std::string& r_new_model_part_name = mThisParameters["new_model_part_name"].GetString();
    ModelPart& r_rigid_model_part = r_new_model_part_name != r_model_part.Name() ? r_model_part.HasSubModelPart(r_new_model_part_name) ? r_model_part.GetSubModelPart(r_new_model_part_name) : r_model_part.CreateSubModelPart(r_new_model_part_name) : r_model_part;

    // Reorder constrains
    IndexType constraint_id = 1;
    for (auto& r_constrain : r_root_model_part.MasterSlaveConstraints()) {
        r_constrain.SetId(constraint_id);
        ++constraint_id;
    }

    // Getting list of variables
    std::vector<Variable<double>> master_double_list_variables, slave_double_list_variables;
    std::vector<VariableComponent<ComponentType>> master_components_list_variables, slave_components_list_variables;
    const std::string& r_master_variable_name = mThisParameters["master_variable_name"].GetString();
    // The master variable
    if(KratosComponents<Variable<double>>::Has(r_master_variable_name)){
        Variable<double> variable = KratosComponents<Variable<double>>::Get(r_master_variable_name);
        master_double_list_variables.push_back(variable);
    } else if (KratosComponents< VariableComponent<ComponentType>>::Has(r_master_variable_name)) {
        VariableComponent<ComponentType> variable = KratosComponents<VariableComponent<ComponentType>>::Get(r_master_variable_name);
        master_components_list_variables.push_back(variable);
    } else if (KratosComponents< Variable< array_1d< double, 3> > >::Has(r_master_variable_name)) {
        VariableComponent<ComponentType> variable_x = KratosComponents<VariableComponent<ComponentType>>::Get(r_master_variable_name+"_X");
        master_components_list_variables.push_back(variable_x);
        VariableComponent<ComponentType> variable_y = KratosComponents<VariableComponent<ComponentType>>::Get(r_master_variable_name+"_Y");
        master_components_list_variables.push_back(variable_y);
        if (r_root_model_part.GetProcessInfo()[DOMAIN_SIZE] == 3) {
            VariableComponent<ComponentType> variable_z = KratosComponents<VariableComponent<ComponentType>>::Get(r_master_variable_name+"_Z");
            master_components_list_variables.push_back(variable_z);
        }
    } else {
        KRATOS_ERROR << "Only double, components and vector variables are allowed in the variables list." ;
    }
    const std::string& r_slave_variable_name = mThisParameters["slave_variable_name"].GetString();
    // We get the slave variable list
    if (r_slave_variable_name != "") {
        if(KratosComponents<Variable<double>>::Has(r_slave_variable_name)){
            Variable<double> variable = KratosComponents<Variable<double>>::Get(r_slave_variable_name);
            slave_double_list_variables.push_back(variable);
        } else if (KratosComponents< VariableComponent<ComponentType>>::Has(r_slave_variable_name)) {
            VariableComponent<ComponentType> variable = KratosComponents<VariableComponent<ComponentType>>::Get(r_slave_variable_name);
            slave_components_list_variables.push_back(variable);
        } else if (KratosComponents< Variable< array_1d< double, 3> > >::Has(r_slave_variable_name)) {
            VariableComponent<ComponentType> variable_x = KratosComponents<VariableComponent<ComponentType>>::Get(r_slave_variable_name+"_X");
            slave_components_list_variables.push_back(variable_x);
            VariableComponent<ComponentType> variable_y = KratosComponents<VariableComponent<ComponentType>>::Get(r_slave_variable_name+"_Y");
            slave_components_list_variables.push_back(variable_y);
            if (r_root_model_part.GetProcessInfo()[DOMAIN_SIZE] == 3) {
                VariableComponent<ComponentType> variable_z = KratosComponents<VariableComponent<ComponentType>>::Get(r_slave_variable_name+"_Z");
                slave_components_list_variables.push_back(variable_z);
            }
        } else {
            KRATOS_ERROR << "Only double, components and vector variables are allowed in the variables list." ;
        }
    } else { // Else we consider exactly the same list of variables
        for (auto& var : master_double_list_variables)
            slave_double_list_variables.push_back(var);
        for (auto& var : master_components_list_variables)
            slave_components_list_variables.push_back(var);
    }

    // Getting index of the master node
    const int master_node_id = mThisParameters["master_node_id"].GetInt();

    // We iterate over the nodes of the rigid model part
    auto nodes_array = r_rigid_model_part.Nodes();
    const int number_of_nodes = static_cast<int>(nodes_array.size());

    // List of variables
    const SizeType number_of_double_variables = master_double_list_variables.size();
    const SizeType number_of_components_variables = master_components_list_variables.size();

    // We get the relation and constant
    const double relation = mThisParameters["relation"].GetDouble();
    const double constant = mThisParameters["constant"].GetDouble();

    // Reference constraint
    const auto& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");

    #pragma omp parallel
    {
        ConstraintContainerType constraints_buffer;
        constraints_buffer.reserve(number_of_nodes);

        // If we master node ID is zero then we get the first node of the model part
        NodeType::Pointer p_master_node = (master_node_id == 0) ? *(r_rigid_model_part.Nodes().begin()).base() : r_root_model_part.pGetNode(master_node_id);
        #pragma omp for
        for (int i = 0; i < number_of_nodes; ++i) {
            auto it_node = nodes_array.begin() + i;
            if (it_node->Id() != p_master_node->Id()) {
                for (IndexType i_var = 0; i_var < number_of_double_variables; ++i_var) {
                    auto constraint = r_clone_constraint.Create(constraint_id + (i * number_of_double_variables + i_var) + 1, *p_master_node, master_double_list_variables[i_var], *it_node, slave_double_list_variables[i_var], relation, constant);
                    (constraints_buffer).insert((constraints_buffer).begin(), constraint);
                }
                for (IndexType i_var = 0; i_var < number_of_components_variables; ++i_var) {
                    auto constraint = r_clone_constraint.Create(constraint_id + (i * number_of_components_variables + i_var) + 1, *p_master_node, master_components_list_variables[i_var], *it_node, slave_components_list_variables[i_var], relation, constant);
                    (constraints_buffer).insert((constraints_buffer).begin(), constraint);
                }
            }
        }

        // We transfer
        #pragma omp critical
        {
            r_rigid_model_part.AddMasterSlaveConstraints(constraints_buffer.begin(),constraints_buffer.end());
            mrThisModelPart.AddMasterSlaveConstraints(constraints_buffer.begin(),constraints_buffer.end());
        }
    }

    KRATOS_CATCH("")
}

// class ImposeRigidMovementProcess
} // namespace Kratos
