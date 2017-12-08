import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest

import math

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return PostProcessEigenvaluesProcess(Model, settings["Parameters"])

class PostProcessEigenvaluesProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):
    def __init__(self, Model, settings):

        default_settings = KratosMultiphysics.Parameters(
            """
            {
                "result_file_name" : "Structure",
                "result_file_format_use_ascii" : false,
                "computing_model_part_name"   : "computing_domain",
                "animation_steps"   :  20,
                "list_of_result_variables" : ["DISPLACEMENT"],
                "label_type" : "frequency"
            }
            """
        );

        settings.ValidateAndAssignDefaults(default_settings)

        KratosMultiphysics.Process.__init__(self)
        model_part = Model[settings["computing_model_part_name"].GetString()]
        settings.RemoveValue("computing_model_part_name")
        self.post_eigen_process = StructuralMechanicsApplication.PostprocessEigenvaluesProcess(
                                    model_part, settings)
                                                                              
    def ExecuteInitialize(self):
        self.post_eigen_process.ExecuteInitialize()
    
    def ExecuteBeforeSolutionLoop(self):
        self.post_eigen_process.ExecuteBeforeSolutionLoop()
    
    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        pass
              
    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        self.post_eigen_process.ExecuteFinalize()
