from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

#Base class to develop other solvers
class ModelManager(object):
    """The base class for solid mechanic model build process.

    This class provides functions for importing and exporting models,
    adding nodal variables and dofs.

    """
    def __init__(self, custom_settings):

        default_settings = KratosMultiphysics.Parameters("""
        {
           "model_name": "solid_domain",
           "dimension": 3,
           "bodies_list": [],
           "domain_parts_list": [],
           "processes_parts_list": [],
           "output_model_part_name": "output_domain",
           "computing_model_part_name": "computing_domain",
           "input_file_settings": {
                "type" : "mdpa",
                "name" : "unknown_name",
                "label": 0
           },
           "variables":[],
           "dofs": []
         }
        """)


        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        self.settings["input_file_settings"].ValidateAndAssignDefaults(default_settings["input_file_settings"])

        # Set void model
        self.main_model_part = self._create_main_model_part()
        self.model = self._create_model()

        # Process Info
        self.process_info = self.main_model_part.ProcessInfo

        # Variables and Dofs settings
        self.nodal_variables = []
        self.dof_variables   = []
        self.dof_reactions   = []


    def ImportModel(self):

        self._add_variables()

        #print("::[Model_Manager]:: Importing model part.")
        problem_path = os.getcwd()
        input_filename = self.settings["input_file_settings"]["name"].GetString()

        if(self.settings["input_file_settings"]["type"].GetString() == "mdpa"):
            # Import model part from mdpa file.
            print("  (reading file: "+ input_filename + ".mdpa)")
            #print("   " + os.path.join(problem_path, input_filename) + ".mdpa ")
            sys.stdout.flush()

            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.SPACE_DIMENSION, self.settings["dimension"].GetInt())
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.settings["dimension"].GetInt()) # Legacy

            KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)

            # Check and prepare computing model part and import constitutive laws.
            self._execute_after_reading()

            self._add_dofs()

            # Somewhere must ask if you want to clean previous files
            self._clean_previous_result_files()

        elif(self.settings["input_file_settings"]["type"].GetString() == "rest"):
            # Import model part from restart file.
            restart_path = os.path.join(problem_path, self.settings["input_file_settings"]["name"].GetString() + "__" + str(self.settings["input_file_settings"]["label"].GetInt() ) )
            if(os.path.exists(restart_path+".rest") == False):
                raise Exception("Restart file not found: " + restart_path + ".rest")
            print("   Loading Restart file: ", restart_path + ".rest ")
            # set serializer flag
            serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_NO_TRACE      # binary
            # serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ERROR # ascii
            # serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ALL   # ascii

            serializer = KratosMultiphysics.Serializer(restart_path, serializer_flag)
            serializer.Load(self.main_model_part.Name, self.main_model_part)

            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True
            #I use it to rebuild the contact conditions.
            load_step = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] +1;
            self.main_model_part.ProcessInfo[KratosMultiphysics.LOAD_RESTART] = load_step
            # print("   Finished loading model part from restart file ")

            computing_model_part = self.settings["computing_model_part_name"].GetString()
            self._add_model_part_to_model(computing_model_part)

            # Get the list of the model_part's in the object Model
            for i in range(self.settings["domain_parts_list"].size()):
                part_name = self.settings["domain_parts_list"][i].GetString()
                self._add_model_part_to_model(part_name)

            for i in range(self.settings["processes_parts_list"].size()):
                part_name = self.settings["processes_parts_list"][i].GetString()
                self._add_model_part_to_model(part_name)

        else:
            raise Exception("Other input options are not yet implemented.")


        #print ("::[Model_Manager]:: Finished importing model part")
        print ("::[Model_Manager]:: Model Ready")


    def ExportModel(self):
        name_out_file = self.settings["input_file_settings"]["name"].GetString()+".out"
        file = open(name_out_file + ".mdpa","w")
        file.close()
        # Model part writing
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

    def CleanModel(self):
        self._clean_body_parts()

    ########

    def GetProcessInfo(self):
        return self.process_info

    def GetModel(self):
        return self.model

    def GetMainModelPart(self):
        return self.main_model_part

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.settings["computing_model_part_name"].GetString())

    def GetOutputModelPart(self):
        #return self.main_model_part.GetSubModelPart(self.settings["output_model_part_name"].GetString())
        return self.main_model_part.GetSubModelPart(self.settings["computing_model_part_name"].GetString())

    def SaveRestart(self):
        pass #one should write the restart file here

    def SetVariables(self, variables):
        self.nodal_variables = self.nodal_variables + variables


    #### Model manager internal methods ####

    def _create_main_model_part(self):
        # Defining the model_part
        main_model_part = KratosMultiphysics.ModelPart(self.settings["model_name"].GetString())
        return main_model_part

    def _create_model(self):
        #TODO: replace this "model" for real one once available in kratos core
        model = {self.settings["model_name"].GetString() : self.main_model_part}
        return model

    def _add_model_part_to_model(self, part_name):
        if( self.main_model_part.HasSubModelPart(part_name) ):
            self.model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})

    def _create_sub_model_part(self, part_name):
        self.main_model_part.CreateSubModelPart(part_name)
        self._add_model_part_to_model(part_name)

    def _add_variables(self):

        self._set_variables()

        self.nodal_variables = self.nodal_variables + self.dof_variables + self.dof_reactions
        self.nodal_variables = list(set(self.nodal_variables))

        self.nodal_variables = [self.nodal_variables[i] for i in range(0,len(self.nodal_variables)) if self.nodal_variables[i] != 'NOT_DEFINED']
        self.nodal_variables.sort()

        for variable in self.nodal_variables:
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.KratosGlobals.GetVariable(variable))
            #print(" Added variable ", KratosMultiphysics.KratosGlobals.GetVariable(variable),"(",variable,")")

        #print(self.nodal_variables)
        #print("::[Model_Manager]:: General Variables ADDED")


    def _add_dofs(self):
        AddDofsProcess = KratosSolid.AddDofsProcess(self.main_model_part, self.dof_variables, self.dof_reactions)
        AddDofsProcess.Execute()

        #print(self.dof_variables + self.dof_reactions)
        #print("::[Model_Manager]:: DOF's ADDED")

    def _set_variables(self):

        # Add input variables supplied
        self._set_input_variables()

        # Add displacements
        self.dof_variables = self.dof_variables + ['DISPLACEMENT']
        self.dof_reactions = self.dof_reactions + ['REACTION']

        # Add dynamic variables
        self.nodal_variables = self.nodal_variables + ['VELOCITY','ACCELERATION']

        # Add specific variables for the problem conditions
        # self.nodal_variables = self.nodal_variables + ['VOLUME_ACCELERATION']
        # they must be added by the process
        # self.nodal_variables = self.nodal_variables + ['POSITIVE_FACE_PRESSURE','NEGATIVE_FACE_PRESSURE','POINT_LOAD','LINE_LOAD','SURFACE_LOAD','POINT_STIFFNESS','LINE_STIFFNESS','SURFACE_STIFFNESS','BALLAST_COEFFICIENT']


        # Add rotational variables
        if self._check_input_dof("ROTATION"):
            # Add specific variables for the problem (rotation dofs)
            self.dof_variables = self.dof_variables + ['ROTATION']
            self.dof_reactions = self.dof_reactions + ['TORQUE']

            self.nodal_variables = self.nodal_variables + ['ANGULAR_VELOCITY','ANGULAR_ACCELERATION']
            # Add large rotation variables
            self.nodal_variables = self.nodal_variables + ['STEP_DISPLACEMENT','STEP_ROTATION','DELTA_ROTATION']
            # Add specific variables for the problem conditions
            # self.nodal_variables = self.nodal_variables + ['POINT_MOMENT','LINE_MOMENT','SURFACE_MOMENT','PLANE_POINT_MOMENT','PLANE_LINE_MOMENT']

        # Add pressure variables
        if self._check_input_dof("PRESSURE"):
            # Add specific variables for the problem (pressure dofs)
            self.dof_variables = self.dof_variables + ['PRESSURE']
            self.dof_reactions = self.dof_reactions + ['PRESSURE_REACTION']

        # Add contat variables
        if self._check_input_dof("LAGRANGE_MULTIPLIER"):
            # Add specific variables for the problem (contact dofs)
            self.dof_variables = self.dof_variables + ['LAGRANGE_MULTIPLIER_NORMAL']
            self.dof_reactions = self.dof_reactions + ['LAGRANGE_MULTIPLIER_NORMAL_REACTION']

        # Add water displacement variables
        if self._check_input_dof("WATER_DISPLACEMENT"):
            self.dof_variables = self.dof_variables + ['WATER_DISPLACEMENT','WATER_VELOCITY','WATER_ACCELERATION']
            self.dof_reactions = self.dof_reactions + ['WATER_DISPLACEMENT_REACTION','WATER_VELOCITY_REACTION','WATER_ACCELERATION_REACTION']

        # Add water pressure variables
        if self._check_input_dof("WATER_PRESSURE"):
            self.dof_variables = self.dof_variables + ['WATER_PRESSURE', 'WATER_PRESSURE_VELOCITY','WATER_PRESSURE_ACCELERATIONN']
            self.dof_reactions = self.dof_reactions + ['REACTION_WATER_PRESSURE', 'WATER_PRESSURE_VELOCITY_REACTION', 'WATER_PRESSURE_ACCELERATION_REACTION']

        # Add jacobian variables
        if self._check_input_dof("JACOBIAN"):
            self.dof_variables = self.dof_variables + ['JACOBIAN']
            self.dof_reactions = self.dof_reactions + ['REACTION_JACOBIAN']


    def _set_input_variables(self):
        variables_list = self.settings["variables"]
        for i in range(0, variables_list.size() ):
            self.nodal_variables.append(variables_list[i].GetString())

    def _check_input_dof(self, variable):
        dofs_list = self.settings["dofs"]
        for i in range(0, dofs_list.size() ):
            if dofs_list[i].GetString() == variable:
                return True
        return False

    #
    def _execute_after_reading(self):

        self._create_sub_model_part(self.settings["computing_model_part_name"].GetString())
        #self._create_sub_model_part(self.settings["output_model_part_name"].GetString())

        # Build bodies
        if( self._has_bodies() ):
            self._build_bodies()

        # Build computing domain
        self._build_computing_domain()

        # Build model
        self._build_model()

    #
    def _build_bodies(self):

        #construct body model parts:
        solid_body_model_parts = []
        fluid_body_model_parts = []
        rigid_body_model_parts = []

        void_flags = KratosSolid.FlagsContainer()

        bodies_list = self.settings["bodies_list"]
        for i in range(bodies_list.size()):
            #create body model part
            body_model_part_name = bodies_list[i]["body_name"].GetString()
            self.main_model_part.CreateSubModelPart(body_model_part_name)
            body_model_part = self.main_model_part.GetSubModelPart(body_model_part_name)

            print("::[Model_Prepare]::Body Created :", body_model_part_name)
            body_model_part.ProcessInfo = self.main_model_part.ProcessInfo
            body_model_part.Properties  = self.main_model_part.Properties

            #build body from their parts
            body_parts_name_list = bodies_list[i]["parts_list"]
            body_parts_list = []
            for j in range(body_parts_name_list.size()):
                body_parts_list.append(self.main_model_part.GetSubModelPart(body_parts_name_list[j].GetString()))

            body_model_part_type = bodies_list[i]["body_type"].GetString()

            for part in body_parts_list:
                entity_type = "Nodes"
                if (body_model_part_type=="Fluid"):
                    assign_flags = KratosSolid.FlagsContainer()
                    assign_flags.PushBack(KratosMultiphysics.FLUID)
                    transfer_process = KratosSolid.TransferEntitiesProcess(body_model_part,part,entity_type,void_flags,assign_flags)
                    transfer_process.Execute()
                elif (body_model_part_type=="Solid"):
                    assign_flags = KratosSolid.FlagsContainer()
                    assign_flags.PushBack(KratosMultiphysics.SOLID)
                    transfer_process = KratosSolid.TransferEntitiesProcess(body_model_part,part,entity_type,void_flags,assign_flags)
                    transfer_process.Execute()
                elif (body_model_part_type=="Rigid"):
                    assign_flags = KratosSolid.FlagsContainer()
                    assign_flags.PushBack(KratosMultiphysics.RIGID)
                    assign_flags.PushBack(KratosMultiphysics.BOUNDARY)
                    transfer_process = KratosSolid.TransferEntitiesProcess(body_model_part,part,entity_type,void_flags,assign_flags)
                    transfer_process.Execute()

                entity_type = "Elements"
                transfer_process = KratosSolid.TransferEntitiesProcess(body_model_part,part,entity_type)
                transfer_process.Execute()
                entity_type = "Conditions"
                transfer_process = KratosSolid.TransferEntitiesProcess(body_model_part,part,entity_type)
                transfer_process.Execute()

            if( body_model_part_type == "Solid" ):
                body_model_part.Set(KratosMultiphysics.SOLID)
                solid_body_model_parts.append(self.main_model_part.GetSubModelPart(body_model_part_name))
            if( body_model_part_type == "Fluid" ):
                body_model_part.Set(KratosMultiphysics.FLUID)
                fluid_body_model_parts.append(self.main_model_part.GetSubModelPart(body_model_part_name))
            if( body_model_part_type == "Rigid" ):
                body_model_part.Set(KratosMultiphysics.RIGID)
                rigid_body_model_parts.append(self.main_model_part.GetSubModelPart(body_model_part_name))

        #add walls in fluid domains:
        transfer_flags = KratosSolid.FlagsContainer()
        transfer_flags.PushBack(KratosMultiphysics.RIGID)
        transfer_flags.PushBack(KratosMultiphysics.NOT_FLUID)

        entity_type = "Nodes"
        for fluid_part in fluid_body_model_parts:
            for rigid_part in rigid_body_model_parts:
                transfer_process = KratosSolid.TransferEntitiesProcess(fluid_part,rigid_part,entity_type,transfer_flags)
                transfer_process.Execute()


    #
    def _build_computing_domain(self):

        # The computing_model_part is labeled 'KratosMultiphysics.ACTIVE' flag (in order to recover it)
        computing_model_part_name  = self.settings["computing_model_part_name"].GetString()
        sub_model_part_names       = self.settings["domain_parts_list"]
        processes_model_part_names = self.settings["processes_parts_list"]

        domain_parts = []
        for i in range(sub_model_part_names.size()):
            domain_parts.append(self.main_model_part.GetSubModelPart(sub_model_part_names[i].GetString()))
        processes_parts = []
        for i in range(processes_model_part_names.size()):
            processes_parts.append(self.main_model_part.GetSubModelPart(processes_model_part_names[i].GetString()))

        computing_model_part = self.main_model_part.GetSubModelPart(computing_model_part_name)
        computing_model_part.ProcessInfo = self.main_model_part.ProcessInfo
        computing_model_part.Properties  = self.main_model_part.Properties

        #set flag to identify the solid model part :: solid application
        computing_model_part.Set(KratosMultiphysics.SOLID)
        #set flag to identify the computing model part
        computing_model_part.Set(KratosMultiphysics.ACTIVE)

        entity_type = "Nodes"
        transfer_process = KratosSolid.TransferEntitiesProcess(computing_model_part,self.main_model_part,entity_type)
        transfer_process.Execute()

        for part in domain_parts:
            entity_type = "Elements"
            transfer_process = KratosSolid.TransferEntitiesProcess(computing_model_part,part,entity_type)
            transfer_process.Execute()

        for part in processes_parts:
            part.Set(KratosMultiphysics.BOUNDARY)
            entity_type = "Conditions"
            #condition flags as BOUNDARY or CONTACT are reserved to composite or contact conditions (do not set it here)
            transfer_process = KratosSolid.TransferEntitiesProcess(computing_model_part,part,entity_type)
            transfer_process.Execute()

    #
    def _build_model(self):

        if( self._has_bodies() ):
            bodies_list = self.settings["bodies_list"]
            for i in range(bodies_list.size()):
                body_parts_name_list = bodies_list[i]["parts_list"]
                for j in range(body_parts_name_list.size()):
                    part_name = body_parts_name_list[j].GetString()
                    if( self.main_model_part.HasSubModelPart(part_name) ):
                        self.model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})
        else:
            # Get the list of the model_part's in the object Model
            for i in range(self.settings["domain_parts_list"].size()):
                part_name = self.settings["domain_parts_list"][i].GetString()
                if( self.main_model_part.HasSubModelPart(part_name) ):
                    self.model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})

        for i in range(self.settings["processes_parts_list"].size()):
            part_name = self.settings["processes_parts_list"][i].GetString()
            if( self.main_model_part.HasSubModelPart(part_name) ):
                self.model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})

    #
    def _clean_body_parts(self):

        #delete body parts: (materials have to be already assigned)
        if( self._has_bodies() ):
            bodies_list = self.settings["bodies_list"]
            for i in range(bodies_list.size()):
                #get body parts
                body_parts_name_list = bodies_list[i]["parts_list"]
                for j in range(body_parts_name_list.size()):
                    self.main_model_part.RemoveSubModelPart(body_parts_name_list[j].GetString())
                    print("::[Model_Prepare]::Body Part Removed:", body_parts_name_list[j].GetString())

    #
    def _has_bodies(self):
        if( self.settings.Has("bodies_list") ):
            if( self.settings["bodies_list"].size() > 0 ):
                return True
        return False

    #
    def _clean_previous_result_files(self):

        file_endings = [".post.bin",".post.msh",".post.res",".post.lst",".post.csv",".rest"]
        problem_path = os.getcwd()
        for file_end in file_endings:
            filelist = [f for f in os.listdir(problem_path) if f.endswith(file_end)]

            for f in filelist:
                try:
                    os.remove(f)
                except OSError:
                    pass
