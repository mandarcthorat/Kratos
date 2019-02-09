from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
from fluid_dynamics_analysis import FluidDynamicsAnalysis
import KratosMultiphysics.kratos_utilities as kratos_utils

try:
    import KratosMultiphysics.ExternalSolversApplication
    have_external_solvers = True
except ImportError as e:
    have_external_solvers = False

import time
import os
import sys

import KratosMultiphysics.KratosUnittest as UnitTest

# Class to navigate through the folders
class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

# Class derived from the UnitTest (KratosMultiphysics.KratosUnittest) class
class TwoFluidInletTest(UnitTest.TestCase):

    def __init__(self):
        self.waterLevel = 0.5
        self.work_folder = "TwoFluidInletTest"
        self.settings = "parameters_serial.json"
        self.check_tolerance = 1e-5
        self.check_toleranceDistance = 0.05
        self.gravitationalAcceleration = 9.81
        self.domainHeight = 1.0
        self.rho1 = 1000.0
        self.rho2 = 1.0
        # switch here for output
        self.print_output = False

    # runs the three dimensional test case
    def run_Serial_Inlet_Test(self):
        with open("TwoFluidInletTest/parameters_serial.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()

            if self.print_output:
                parameters["output_processes"].AddValue("gid_output", KratosMultiphysics.Parameters(R'''[{
                    "python_module" : "gid_output_process",
                    "kratos_module" : "KratosMultiphysics",
                    "process_name"  : "GiDOutputProcess",
                    "help"          : "This process writes postprocessing files for GiD",
                    "Parameters"    : {
                        "model_part_name"        : "FluidModelPart",
                        "output_name"            : "FluidModelPart",
                        "postprocess_parameters" : {
                            "result_file_configuration" : {
                                "gidpost_flags"       : {
                                    "GiDPostMode"           : "GiD_PostBinary",
                                    "WriteDeformedMeshFlag" : "WriteDeformed",
                                    "WriteConditionsFlag"   : "WriteConditions",
                                    "MultiFileFlag"         : "SingleFile"
                                },
                                "file_label"          : "time",
                                "output_control_type" : "time",
                                "output_frequency"    : 0.1,
                                "body_output"         : true,
                                "node_output"         : false,
                                "skin_output"         : false,
                                "plane_output"        : [],
                                "nodal_results"       : ["VELOCITY","PRESSURE","DISTANCE","DENSITY","DYNAMIC_VISCOSITY"],
                                "gauss_point_results" : []
                            },
                            "point_data_configuration"  : []
                        }
                    }
                }]''') )

            # running
            self.simulation = FluidDynamicsAnalysisWithFlush3D(model,parameters)
            self.simulation.Run()
            self._check_results_serial()
            self._clean_up()


    # runs the three dimensional test case
    def run_MPI_Inlet_Test(self):
        with open("TwoFluidInletTest/parameters_mpi.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()

            if self.print_output:
                parameters["output_processes"].AddValue("gid_output", KratosMultiphysics.Parameters(R'''[{
                    "python_module" : "gid_output_process",
                    "kratos_module" : "KratosMultiphysics",
                    "process_name"  : "GiDOutputProcess",
                    "help"          : "This process writes postprocessing files for GiD",
                    "Parameters"    : {
                        "model_part_name"        : "FluidModelPart",
                        "output_name"            : "FluidModelPart",
                        "postprocess_parameters" : {
                            "result_file_configuration" : {
                                "gidpost_flags"       : {
                                    "GiDPostMode"           : "GiD_PostBinary",
                                    "WriteDeformedMeshFlag" : "WriteDeformed",
                                    "WriteConditionsFlag"   : "WriteConditions",
                                    "MultiFileFlag"         : "SingleFile"
                                },
                                "file_label"          : "time",
                                "output_control_type" : "time",
                                "output_frequency"    : 0.1,
                                "body_output"         : true,
                                "node_output"         : false,
                                "skin_output"         : false,
                                "plane_output"        : [],
                                "nodal_results"       : ["VELOCITY","PRESSURE","DISTANCE","DENSITY","DYNAMIC_VISCOSITY"],
                                "gauss_point_results" : []
                            },
                            "point_data_configuration"  : []
                        }
                    }
                }]''') )

            # running
            self.simulation = FluidDynamicsAnalysisWithFlush3D(model,parameters)
            self.simulation.Run()
            self._check_results_mpi()
            self._clean_up()



    def _check_results_serial( self ):

        model = self.simulation._GetSolver().GetComputingModelPart().GetModel()
        inlet = model.GetModelPart("AutomaticInlet3D_Inlet")

        for node in inlet.Nodes:

            # checking the formula-based velocity field at the inlet
            velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            if (node.Z < 0.7):
                v_x_ref = node.Z * 2.0
                self.assertAlmostEqual( v_x_ref, velocity[0], delta = self.check_tolerance)
            if (node.Z > 0.7):
                v_x_ref = 0.1
                self.assertAlmostEqual( v_x_ref, velocity[0], delta = self.check_tolerance)
            self.assertAlmostEqual(0.0, velocity[1], delta = self.check_tolerance)
            self.assertAlmostEqual(0.0, velocity[2], delta = self.check_tolerance)

        test_node = self.simulation._GetSolver().GetComputingModelPart().GetNode( 6 )
        self.assertAlmostEqual( -0.125, test_node.GetSolutionStepValue(KratosMultiphysics.DISTANCE), delta = self.check_toleranceDistance)
        test_node = self.simulation._GetSolver().GetComputingModelPart().GetNode( 12 )
        self.assertAlmostEqual( -0.125, test_node.GetSolutionStepValue(KratosMultiphysics.DISTANCE), delta = self.check_toleranceDistance)
        test_node = self.simulation._GetSolver().GetComputingModelPart().GetNode( 43 )
        self.assertAlmostEqual( -0.3, test_node.GetSolutionStepValue(KratosMultiphysics.DISTANCE), delta = self.check_toleranceDistance)
        test_node = self.simulation._GetSolver().GetComputingModelPart().GetNode( 100 )
        self.assertAlmostEqual( 0.28749998375, test_node.GetSolutionStepValue(KratosMultiphysics.DISTANCE), delta = self.check_toleranceDistance)



    def _check_results_mpi( self ):

        model = self.simulation._GetSolver().GetComputingModelPart().GetModel()
        inlet = model.GetModelPart("AutomaticInlet3D_Inlet")

        for node in inlet.Nodes:

            # checking the formula-based velocity field at the inlet
            velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            if (node.Z < 0.7):
                v_x_ref = node.Z * 2.0
                self.assertAlmostEqual( v_x_ref, velocity[0], delta = self.check_tolerance)
            if (node.Z > 0.7):
                v_x_ref = 0.1
                self.assertAlmostEqual( v_x_ref, velocity[0], delta = self.check_tolerance)
            self.assertAlmostEqual(0.0, velocity[1], delta = self.check_tolerance)
            self.assertAlmostEqual(0.0, velocity[2], delta = self.check_tolerance)


    def _clean_up( self ):

        kratos_utils.DeleteFileIfExisting("FluidModelPart.post.bin")
        kratos_utils.DeleteFileIfExisting("tests.post.lst")
        kratos_utils.DeleteFileIfExisting("tests.post.bin")

        for i in range(0,100):
            file_name = "TwoFluidInletTest/test_inlet_" + str(i) + ".mdpa"
            kratos_utils.DeleteFileIfExisting(file_name)

        print( "clean" )


class FluidDynamicsAnalysisWithFlush3D(FluidDynamicsAnalysis):

    def __init__(self,model,project_parameters,flush_frequency=10.0):
        super(FluidDynamicsAnalysisWithFlush3D,self).__init__(model,project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def ModifyInitialGeometry(self):
        init_h = 0.3
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            distance = node.Z - init_h
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)

    def FinalizeSolutionStep(self):
        super(FluidDynamicsAnalysisWithFlush3D,self).FinalizeSolutionStep()
        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now


if __name__ == "__main__":

    test = TwoFluidInletTest()

    test.run_Serial_Inlet_Test()

    # test.run_MPI_Inlet_Test()