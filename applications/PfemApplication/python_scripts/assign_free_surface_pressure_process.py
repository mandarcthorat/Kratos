from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemApplication as KratosPfem

def Factory(custom_settings, Model):
    if( not isinstance(custom_settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignFreeSurfacePressureProcess(Model, custom_settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignFreeSurfacePressureProcess(KratosMultiphysics.Process):
    def __init__(self, Model, custom_settings ):
        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
             "model_part_name": "main_domain"
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model = Model

    @classmethod
    def GetVariables(self):
        nodal_variables = []
        return nodal_variables

    def ExecuteInitialize(self):

        # set model part
        self.model_part = self.model[self.settings["model_part_name"].GetString()]

    def ExecuteInitializeSolutionStep(self):
        for i in self.model_part.Nodes:
            i.Free(KratosMultiphysics.PRESSURE)
            # if i.Is(KratosMultiphysics.FREE_SURFACE) or (i.Is(KratosMultiphysics.RIGID) and i.IsNot(KratosMultiphysics.FLUID)):
            if (i.Is(KratosMultiphysics.RIGID) and i.IsNot(KratosMultiphysics.FLUID)):
                i.Fix(KratosMultiphysics.PRESSURE)
                i.SetSolutionStepValue(KratosMultiphysics.PRESSURE,0.0)
                i.SetSolutionStepValue(KratosMultiphysics.PRESSURE,1,0.0)
