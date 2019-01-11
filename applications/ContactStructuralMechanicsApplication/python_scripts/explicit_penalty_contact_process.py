from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KM

import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

def Factory(settings, Model):
    if(type(settings) != KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ExplicitPenaltyContactProcess(Model, settings["Parameters"])

import sys

import penalty_contact_process

class ExplicitPenaltyContactProcess(penalty_contact_process.PenaltyContactProcess):
    """This class is used in order to compute the contact using a mortar penalty formulation in explicit integration schemes

    This class constructs the model parts containing the contact conditions and
    initializes parameters and variables related with the contact. The class creates
    search utilities to be used to create the contact pairs

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """

        # NOTE: Due to recursive check "contact_model_part" and "assume_master_slave" requires to pre-define configurations, if more that 10 pairs of contact are required, just add. I assume nobody needs that much
        # Settings string in json format
        default_parameters = KM.Parameters("""
        {
            "help"                        : "This class is used in order to compute the contact using a mortar ALM formulation. This class constructs the model parts containing the contact conditions and initializes parameters and variables related with the contact. The class creates search utilities to be used to create the contact pairs",
            "mesh_id"                     : 0,
            "model_part_name"             : "Structure",
            "computing_model_part_name"   : "computing_domain",
            "contact_model_part"          : {"0":[],"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[]},
            "assume_master_slave"         : {"0":[],"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[]},
            "contact_property_ids"        : {"0": 0,"1": 0,"2": 0,"3": 0,"4": 0,"5": 0,"6": 0,"7": 0,"8": 0,"9": 0},
            "contact_type"                : "Frictionless",
            "interval"                    : [0.0,"End"],
            "normal_variation"            : "no_derivatives_computation",
            "frictional_law"              : "Coulomb",
            "tangent_factor"              : 1.0e-1,
            "integration_order"           : 2,
            "clear_inactive_for_post"     : true,
            "search_parameters" : {
                "type_search"                         : "in_radius",
                "simple_search"                       : false,
                "adapt_search"                        : false,
                "search_factor"                       : 3.5,
                "active_check_factor"                 : 0.01,
                "max_number_results"                  : 1000,
                "bucket_size"                         : 4,
                "dynamic_search"                      : false,
                "static_check_movement"               : false,
                "database_step_update"                : 1,
                "consider_gap_threshold"              : false,
                "debug_mode"                          : false,
                "check_gap"                           : "check_mapping"
            },
            "advance_ALM_parameters" : {
                "manual_ALM"                  : false,
                "stiffness_factor"            : 1.0,
                "penalty"                     : 1.0e16,
                "adapt_penalty"               : false,
                "max_gap_factor"              : 1.0e-3
            },
            "alternative_formulations" : {
                "axisymmetric"                : false
            }
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.contact_settings = settings
        self.contact_settings.RecursivelyValidateAndAssignDefaults(default_parameters)

        # Construct the base process.
        super(ExplicitPenaltyContactProcess, self).__init__(Model, self.contact_settings)

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(ExplicitPenaltyContactProcess, self).ExecuteInitialize()

        # Create the dynamic factor process
        self.dynamic_factor_process = CSMA.ComputeDynamicFactorProcess(self.main_model_part)

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super(ExplicitPenaltyContactProcess, self).ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(ExplicitPenaltyContactProcess, self).ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super(ExplicitPenaltyContactProcess, self).ExecuteFinalizeSolutionStep()

        # Calling for the active set utilities (to activate deactivate nodes)
        if self.contact_settings["contact_type"].GetString() == "Frictionless":
            CSMA.ActiveSetUtilities.ComputePenaltyFrictionlessActiveSet(self.computing_model_part)
        else:
            CSMA.ActiveSetUtilities.ComputePenaltyFrictionalActiveSet(self.computing_model_part)

        # Updating value of weighted gap
        KM.VariableUtils().SetVariable(CSMA.WEIGHTED_GAP, 0.0, self.computing_model_part.Nodes);
        CSMA.ContactUtilities.ComputeExplicitContributionConditions(self.computing_model_part)

        # Update the dynamic factors
        self.dynamic_factor_process.Execute()

    def ExecuteBeforeOutputStep(self):
        """ This method is executed right before the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(ExplicitPenaltyContactProcess, self).ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        """ This method is executed right after the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(ExplicitPenaltyContactProcess, self).ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        """ This method is executed in order to finalize the current computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(ExplicitPenaltyContactProcess, self).ExecuteFinalize()

    def _compute_search(self):
        """ This method return if the serach must be computed

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # TODO: Adding a proper check for explicit computations
        if self._get_if_is_interval():
            self.database_step += 1
            global_step = self.main_model_part.ProcessInfo[KM.STEP]
            database_step_update = self.settings["search_parameters"]["database_step_update"].GetInt()
            if self.database_step >= database_step_update or global_step == 1:
                return True
            else:
                return False

    def _initialize_problem_parameters(self):
        """ This method initializes the ALM parameters from the process info

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process (in fact not, to avoid writing twice the values)
        #super(PenaltyContactProcess, self)._initialize_problem_parameters()

        # We call the process info
        process_info = self.main_model_part.ProcessInfo

        if not self.contact_settings["advance_ALM_parameters"]["manual_ALM"].GetBool():
            # We compute NODAL_H that can be used in the search and some values computation
            self.find_nodal_h = KM.FindNodalHProcess(self.computing_model_part)
            self.find_nodal_h.Execute()

            # Computing the scale factors or the penalty parameters (StiffnessFactor * E_mean/h_mean)
            alm_var_parameters = KM.Parameters("""{}""")
            alm_var_parameters.AddValue("stiffness_factor", self.contact_settings["advance_ALM_parameters"]["stiffness_factor"])
            alm_var_parameters.AddValue("penalty_scale_factor", self.contact_settings["advance_ALM_parameters"]["penalty_scale_factor"])
            self.alm_var_process = CSMA.ALMVariablesCalculationProcess(self._get_process_model_part(), KM.NODAL_H, alm_var_parameters)
            self.alm_var_process.Execute()
            process_info[KM.INITIAL_PENALTY] = 1.0e0 * process_info[KM.INITIAL_PENALTY] # We rescale, the process is designed for ALM formulation
        else:
            # We set the values in the process info
            process_info[KM.INITIAL_PENALTY] = self.contact_settings["advance_ALM_parameters"]["penalty"].GetDouble()

        # We set a minimum value
        if process_info[KM.INITIAL_PENALTY] < sys.float_info.epsilon:
            process_info[KM.INITIAL_PENALTY] = 1.0e12

        # We print the parameters considered
        KM.Logger.PrintInfo("INITIAL_PENALTY: ", "{:.2e}".format(process_info[KM.INITIAL_PENALTY]))
