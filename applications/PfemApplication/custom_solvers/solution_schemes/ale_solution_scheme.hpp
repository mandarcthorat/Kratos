//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:        IPouplana $
//   Last modified by:    $Co-Author:   JMCarbonell $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_ALE_SOLUTION_SCHEME_H_INCLUDED )
#define  KRATOS_ALE_SOLUTION_SCHEME_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"
#include "includes/node.h"
#include "geometries/geometry.h"

#include "pfem_application_variables.h"

#include "custom_solvers/solution_schemes/dynamic_scheme.hpp"
#include "custom_utilities/slip_coordinate_transformation.hpp"


namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Bossak time scheme for the incompressible flow problem.
/** This class provides a second order time scheme of the generalized-alpha Newmark
    family of methods. It also includes code required to implement slip conditions
    on the incompressible flow problem and provides the possibility of using a RANS
    model by passing a turbulence model as an argument to the constructor.
    This time scheme is intended to be used in combination with elements of type
    ASGS2D, ASGS3D, VMS or derived classes.
    To use the slip condition SLIP flag must be assigned to nodes. To use
    a wall law in combination with the slip condition, use MonolithicWallCondition to
    mesh the boundary
    @see ASGS2D, ASGS3D, VMS, MonolithicWallConditon
*/

template<class TSparseSpace, class TDenseSpace>
class AleSolutionScheme : public DynamicScheme<TSparseSpace, TDenseSpace>
{
 public:

  ///@name Type Definitions
  ///@{
  KRATOS_CLASS_POINTER_DEFINITION(AleSolutionScheme);

  typedef SolutionScheme<TSparseSpace,TDenseSpace>                             BaseType;
  typedef typename BaseType::SolutionSchemePointerType                  BasePointerType;
  typedef typename BaseType::LocalFlagType                                LocalFlagType;

  typedef DynamicScheme<TSparseSpace,TDenseSpace>                           DerivedType;

  typedef typename BaseType::NodeType                                          NodeType;
  typedef typename BaseType::DofsArrayType                                DofsArrayType;
  typedef typename BaseType::SystemMatrixType                          SystemMatrixType;
  typedef typename BaseType::SystemVectorType                          SystemVectorType;
  typedef typename BaseType::LocalSystemVectorType                LocalSystemVectorType;
  typedef typename BaseType::LocalSystemMatrixType                LocalSystemMatrixType;

  typedef typename BaseType::NodesContainerType                      NodesContainerType;
  typedef typename BaseType::ElementsContainerType                ElementsContainerType;
  typedef typename BaseType::ConditionsContainerType            ConditionsContainerType;

  typedef typename BaseType::IntegrationMethodsVectorType  IntegrationMethodsVectorType;
  typedef typename BaseType::IntegrationMethodsScalarType  IntegrationMethodsScalarType;

  typedef Geometry<NodeType>      GeometryType;
  typedef typename GeometryType::SizeType      SizeType;


  ///@}
  ///@name Life Cycle
  ///@{

    /// Constructor.
    AleSolutionScheme(IntegrationMethodsVectorType& rTimeVectorIntegrationMethods,
                  IntegrationMethodsScalarType& rTimeScalarIntegrationMethods)
        :DerivedType(rTimeVectorIntegrationMethods, rTimeScalarIntegrationMethods)
        , mRotationTool()
    {
    }

    /// Constructor.
    AleSolutionScheme(IntegrationMethodsVectorType& rTimeVectorIntegrationMethods,
                  IntegrationMethodsScalarType& rTimeScalarIntegrationMethods,
                  Flags& rOptions)
        :DerivedType(rTimeVectorIntegrationMethods, rTimeScalarIntegrationMethods, rOptions)
        , mRotationTool()
    {
    }

  // Class SolutionScheme
  ~AleSolutionScheme() override {}

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  /**
   * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
   * @details This is intended to be called just once when the strategy is initialized
   */
  void Initialize(ModelPart& rModelPart) override
  {
    KRATOS_TRY

    DerivedType::Initialize(rModelPart);

    const SizeType dimension = rModelPart.GetProcessInfo()[SPACE_DIMENSION];

    // Set dofs size
    unsigned int DofsSize = 0;
    for(typename IntegrationMethodsVectorType::iterator it=(this->mTimeVectorIntegrationMethods).begin();
        it!=(this->mTimeVectorIntegrationMethods).end(); ++it)
      DofsSize += dimension;

    for(typename IntegrationMethodsScalarType::iterator it=(this->mTimeScalarIntegrationMethods).begin();
        it!=(this->mTimeScalarIntegrationMethods).end(); ++it)
      ++DofsSize;

    mRotationTool.SetBlockSize(DofsSize);

    // Set time integration of MESH_ variables needed for ALE
    for(typename IntegrationMethodsVectorType::iterator it=(this->mTimeVectorIntegrationMethods).begin(); it!=(this->mTimeVectorIntegrationMethods).end(); ++it)
    {
      if( (*it)->GetPrimaryVariableName() == "VELOCITY" )
      {
        (this->mTimeVectorIntegrationMethods).push_back( (*it)->Clone() );
        (this->mTimeVectorIntegrationMethods).back()->SetVariables(MESH_DISPLACEMENT,MESH_VELOCITY,MESH_ACCELERATION,MESH_VELOCITY);
        break;
      }
    }

    KRATOS_CATCH("")
  }

  /**
     Performing the update of the solution.
  */
  //***************************************************************************

  /**
   * Performing the update of the solution
   * @param rModelPart: The model of the problem to solve
   * @param rDofSet: Set of all primary variables
   * @param rDx: incremental update of primary variables
   */
  void Update(ModelPart& rModelPart,
              DofsArrayType& rDofSet,
              SystemVectorType& rDx) override
  {
    KRATOS_TRY;

    mRotationTool.RotateVelocities(rModelPart);

    this->UpdateDofs(rModelPart,rDofSet,rDx);

    // TODO: We must rotate and recover mesh_velocity too (modify mRotationTool methods)
    // if( this->mOptions.Is(LocalFlagType::MOVE_MESH) )
      this->UpdateMeshVelocity(rModelPart,rDx);

    mRotationTool.RecoverVelocities(rModelPart);

    this->UpdateVariables(rModelPart);

    this->MoveMesh(rModelPart,"MESH_DISPLACEMENT");

    KRATOS_CATCH( "" );
  }


  /**
   * @brief Performing the update of the mesh velocity except for slip conditions
   * for slip conditions, mesh velocity must be imposed through a process
   * @details this function must be called only once per iteration
   */
  void UpdateMeshVelocity(ModelPart& rModelPart,
                          SystemVectorType& rDx)
  {
    KRATOS_TRY

    // Updating time derivatives (nodally for efficiency)
    const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector NodePartition;
    OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

    const int nnodes = static_cast<int>(rModelPart.Nodes().size());
    typename NodesContainerType::iterator NodeBegin = rModelPart.Nodes().begin();

#pragma omp parallel for firstprivate(NodeBegin)
    for(int i = 0;  i < nnodes; i++)
    {
      typename NodesContainerType::iterator itNode = NodeBegin + i;

      if(itNode->IsNot(SLIP))
      {
        if (itNode->GetDof(VELOCITY_X).IsFree() )
        {
          itNode->FastGetSolutionStepValue(MESH_VELOCITY_X) += TSparseSpace::GetValue(rDx,itNode->GetDof(VELOCITY_X).EquationId());
        }
        else
        {
          itNode->FastGetSolutionStepValue(MESH_VELOCITY_X) = itNode->FastGetSolutionStepValue(VELOCITY_X);
        }
        if (itNode->GetDof(VELOCITY_Y).IsFree() )
        {
          itNode->FastGetSolutionStepValue(MESH_VELOCITY_Y) += TSparseSpace::GetValue(rDx,itNode->GetDof(VELOCITY_Y).EquationId());
        }
        else
        {
          itNode->FastGetSolutionStepValue(MESH_VELOCITY_Y) = itNode->FastGetSolutionStepValue(VELOCITY_Y);
        }
        if(rModelPart.GetProcessInfo()[SPACE_DIMENSION]==3)
        {
          if (itNode->GetDof(VELOCITY_Z).IsFree() )
          {
            itNode->FastGetSolutionStepValue(MESH_VELOCITY_Z) += TSparseSpace::GetValue(rDx,itNode->GetDof(VELOCITY_Z).EquationId());
          }
          else
          {
            itNode->FastGetSolutionStepValue(MESH_VELOCITY_Z) = itNode->FastGetSolutionStepValue(VELOCITY_Z);
          }
        }
      }
    }

    KRATOS_CATCH("")
  }

  //***************************************************************************

  /**
   * Performing the prediction of the solution
   * @param rModelPart: The model of the problem to solve
   * @param rDofSet set of all primary variables
   * @param rDx: Incremental update of primary variables
   */

  void Predict(ModelPart& rModelPart,
               DofsArrayType& rDofSet,
               SystemVectorType& rDx) override
  {
    KRATOS_TRY;

    this->PredictVariables(rModelPart);

    this->MoveMesh(rModelPart);

    KRATOS_CATCH( "" );
  }

  //***************************************************************************

  /** this function is designed to be called in the builder and solver
      to introduce
      the selected time integration scheme. It "asks" the matrix needed to
      the element and
      performs the operations needed to introduce the seected time
      integration scheme.

      this function calculates at the same time the contribution to the
      LHS and to the RHS
      of the system
  */
  void CalculateSystemContributions(Element::Pointer pCurrentElement,
                                    LocalSystemMatrixType& rLHS_Contribution,
                                    LocalSystemVectorType& rRHS_Contribution,
                                    Element::EquationIdVectorType& EquationId,
                                    ProcessInfo& rCurrentProcessInfo) override
  {
    KRATOS_TRY

    int thread = OpenMPUtils::ThisThread();

    //Initializing the non linear iteration for the current element

    (pCurrentElement) -> CalculateLocalSystem(rLHS_Contribution,rRHS_Contribution, rCurrentProcessInfo);

    (pCurrentElement) -> EquationIdVector(EquationId, rCurrentProcessInfo);

    if ( rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true ){
      (pCurrentElement) -> CalculateSecondDerivativesContributions(this->mMatrix.M[thread],this->mVector.a[thread],rCurrentProcessInfo);
      (pCurrentElement) -> CalculateFirstDerivativesContributions(this->mMatrix.D[thread],this->mVector.v[thread],rCurrentProcessInfo);

      this->AddDynamicTangentsToLHS(rLHS_Contribution,this->mMatrix.D[thread],this->mMatrix.M[thread],rCurrentProcessInfo);
      this->AddDynamicForcesToRHS(rRHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);
    }
    else{
      (pCurrentElement) -> CalculateMassMatrix(this->mMatrix.M[thread], rCurrentProcessInfo);
      (pCurrentElement) -> CalculateLocalVelocityContribution(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

      this->AddDynamicsToLHS (rLHS_Contribution, this->mMatrix.D[thread], this->mMatrix.M[thread], rCurrentProcessInfo);
      this->AddDynamicsToRHS (pCurrentElement, rRHS_Contribution, this->mMatrix.D[thread], this->mMatrix.M[thread], rCurrentProcessInfo);

    }

    // If there is a slip condition, apply it on a rotated system of coordinates
    mRotationTool.Rotate(rLHS_Contribution,rRHS_Contribution,pCurrentElement->GetGeometry());
    mRotationTool.ApplySlipCondition(rLHS_Contribution,rRHS_Contribution,pCurrentElement->GetGeometry());

    KRATOS_CATCH("")

  }

  void Calculate_RHS_Contribution(Element::Pointer pCurrentElement,
                                  LocalSystemVectorType& rRHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  ProcessInfo& rCurrentProcessInfo) override
  {
    KRATOS_TRY

    int thread = OpenMPUtils::ThisThread();

    //Initializing the non linear iteration for the current element

    // Basic operations for the element considered
    (pCurrentElement) -> CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);

    (pCurrentElement) -> EquationIdVector(EquationId, rCurrentProcessInfo);

    if ( rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true ){

      (pCurrentElement) -> CalculateSecondDerivativesRHS(this->mVector.a[thread],rCurrentProcessInfo);
      (pCurrentElement) -> CalculateFirstDerivativesRHS(this->mVector.v[thread],rCurrentProcessInfo);

      this->AddDynamicForcesToRHS(rRHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);

    }
    else{

      (pCurrentElement) -> CalculateMassMatrix(this->mMatrix.M[thread], rCurrentProcessInfo);
      (pCurrentElement) -> CalculateLocalVelocityContribution(this->mMatrix.D[thread], rRHS_Contribution, rCurrentProcessInfo);

      this->AddDynamicsToRHS (pCurrentElement, rRHS_Contribution, this->mMatrix.D[thread], this->mMatrix.M[thread], rCurrentProcessInfo);
    }

    // If there is a slip condition, apply it on a rotated system of coordinates
    mRotationTool.Rotate(rRHS_Contribution,pCurrentElement->GetGeometry());
    mRotationTool.ApplySlipCondition(rRHS_Contribution,pCurrentElement->GetGeometry());

    KRATOS_CATCH("")
  }

  /** functions totally analogous to the precedent but applied to
      the "condition" objects
  */
  void Condition_CalculateSystemContributions(Condition::Pointer pCurrentCondition,
                                              LocalSystemMatrixType& rLHS_Contribution,
                                              LocalSystemVectorType& rRHS_Contribution,
                                              Element::EquationIdVectorType& EquationId,
                                              ProcessInfo& rCurrentProcessInfo) override
  {
    KRATOS_TRY

    int thread = OpenMPUtils::ThisThread();

    // Basic operations for the condition considered
    (pCurrentCondition) -> CalculateLocalSystem(rLHS_Contribution,rRHS_Contribution, rCurrentProcessInfo);

    (pCurrentCondition) -> EquationIdVector(EquationId, rCurrentProcessInfo);

    if ( rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true ){
      (pCurrentCondition) -> CalculateSecondDerivativesContributions(this->mMatrix.M[thread],this->mVector.a[thread],rCurrentProcessInfo);
      (pCurrentCondition) -> CalculateFirstDerivativesContributions(this->mMatrix.D[thread],this->mVector.v[thread],rCurrentProcessInfo);

      this->AddDynamicTangentsToLHS(rLHS_Contribution,this->mMatrix.D[thread],this->mMatrix.M[thread],rCurrentProcessInfo);
      this->AddDynamicForcesToRHS(rRHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);
    }
    else{
      (pCurrentCondition) -> CalculateMassMatrix(this->mMatrix.M[thread], rCurrentProcessInfo);
      (pCurrentCondition) -> CalculateLocalVelocityContribution(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

      this->AddDynamicsToLHS  (rLHS_Contribution, this->mMatrix.D[thread], this->mMatrix.M[thread], rCurrentProcessInfo);
      this->AddDynamicsToRHS  (pCurrentCondition, rRHS_Contribution, this->mMatrix.D[thread], this->mMatrix.M[thread], rCurrentProcessInfo);
    }

    // Rotate contributions (to match coordinates for slip conditions)
    mRotationTool.Rotate(rLHS_Contribution,rRHS_Contribution,pCurrentCondition->GetGeometry());
    mRotationTool.ApplySlipCondition(rLHS_Contribution,rRHS_Contribution,pCurrentCondition->GetGeometry());

    KRATOS_CATCH("")
  }

  void Condition_Calculate_RHS_Contribution(Condition::Pointer pCurrentCondition,
                                            LocalSystemVectorType& rRHS_Contribution,
                                            Element::EquationIdVectorType& EquationId,
                                            ProcessInfo& rCurrentProcessInfo) override
  {
    KRATOS_TRY

    int thread = OpenMPUtils::ThisThread();

    // Basic operations for the condition considered
    (pCurrentCondition) -> CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);

    (pCurrentCondition) -> EquationIdVector(EquationId, rCurrentProcessInfo);

    if ( rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true ){
      (pCurrentCondition) -> CalculateSecondDerivativesRHS(this->mVector.a[thread],rCurrentProcessInfo);
      (pCurrentCondition) -> CalculateFirstDerivativesRHS(this->mVector.v[thread],rCurrentProcessInfo);

      this->AddDynamicForcesToRHS(rRHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);
    }
    else{
      (pCurrentCondition) -> CalculateMassMatrix(this->mMatrix.M[thread], rCurrentProcessInfo);
      (pCurrentCondition) -> CalculateLocalVelocityContribution(this->mMatrix.D[thread], rRHS_Contribution,rCurrentProcessInfo);

      this->AddDynamicsToRHS  (pCurrentCondition, rRHS_Contribution, this->mMatrix.D[thread], this->mMatrix.M[thread], rCurrentProcessInfo);
    }

    // Rotate contributions (to match coordinates for slip conditions)
    mRotationTool.Rotate(rRHS_Contribution,pCurrentCondition->GetGeometry());
    mRotationTool.ApplySlipCondition(rRHS_Contribution,pCurrentCondition->GetGeometry());

    KRATOS_CATCH("")
  }

  // TODO

  // void InitializeNonLinIteration(ModelPart& r_model_part,
  //                                        TSystemMatrixType& A,
  //                                        TSystemVectorType& Dx,
  //                                        TSystemVectorType& b) override
  // {
  //     KRATOS_TRY

  //     if (mpTurbulenceModel != 0) // If not null
  //         mpTurbulenceModel->Execute();

  //     KRATOS_CATCH("")
  // }

  // void FinalizeNonLinIteration(ModelPart &rModelPart, TSystemMatrixType &A, TSystemVectorType &Dx, TSystemVectorType &b) override
  // {
  //     ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

  //     //if orthogonal subscales are computed
  //     if (CurrentProcessInfo[OSS_SWITCH] == 1.0) {
  //         if (rModelPart.GetCommunicator().MyPID() == 0)
  //             std::cout << "Computing OSS projections" << std::endl;


  //         const int nnodes = static_cast<int>(rModelPart.Nodes().size());
  //         auto nbegin = rModelPart.NodesBegin();
  //         #pragma omp parallel for firstprivate(nbegin,nnodes)
  //         for(int i=0; i<nnodes; ++i)
  //         {
  //             auto ind = nbegin + i;
  //             noalias(ind->FastGetSolutionStepValue(ADVPROJ)) = ZeroVector(3);

  //             ind->FastGetSolutionStepValue(DIVPROJ) = 0.0;

  //             ind->FastGetSolutionStepValue(NODAL_AREA) = 0.0;


  //         }//end of loop over nodes

  //         //loop on nodes to compute ADVPROJ   CONVPROJ NODALAREA
  //         array_1d<double, 3 > output = ZeroVector(3);

  //         const int nel = static_cast<int>(rModelPart.Elements().size());
  //         auto elbegin = rModelPart.ElementsBegin();
  //         #pragma omp parallel for firstprivate(elbegin,nel,output)
  //         for(int i=0; i<nel; ++i)
  //         {
  //             auto elem = elbegin + i;
  //             elem->Calculate(ADVPROJ, output, CurrentProcessInfo);
  //         }

  //         rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);
  //         rModelPart.GetCommunicator().AssembleCurrentData(DIVPROJ);
  //         rModelPart.GetCommunicator().AssembleCurrentData(ADVPROJ);

  //         // Correction for periodic conditions
  //         this->PeriodicConditionProjectionCorrection(rModelPart);

  //         #pragma omp parallel for firstprivate(nbegin,nnodes)
  //         for(int i=0; i<nnodes; ++i)
  //         {
  //             auto ind = nbegin + i;
  //             if (ind->FastGetSolutionStepValue(NODAL_AREA) == 0.0)
  //             {
  //                 ind->FastGetSolutionStepValue(NODAL_AREA) = 1.0;
  //                 //KRATOS_WATCH("*********ATTENTION: NODAL AREA IS ZERRROOOO************");
  //             }
  //             const double Area = ind->FastGetSolutionStepValue(NODAL_AREA);
  //             ind->FastGetSolutionStepValue(ADVPROJ) /= Area;
  //             ind->FastGetSolutionStepValue(DIVPROJ) /= Area;
  //         }
  //     }
  // }

  // void FinalizeSolutionStep(ModelPart &rModelPart, TSystemMatrixType &A, TSystemVectorType &Dx, TSystemVectorType &b) override
  // {
  //     Element::EquationIdVectorType EquationId;
  //     LocalSystemVectorType RHS_Contribution;
  //     LocalSystemMatrixType LHS_Contribution;
  //     ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

  //     //for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); ++itNode)

  //     #pragma omp parallel for
  //     for(int k = 0; k<static_cast<int>(rModelPart.Nodes().size()); k++)
  //     {
  //         auto itNode = rModelPart.NodesBegin() + k;
  //         (itNode->FastGetSolutionStepValue(REACTION)).clear();

  //         // calculating relaxed acceleration
  //         const array_1d<double, 3 > & CurrentAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 0);
  //         const array_1d<double, 3 > & OldAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 1);
  //         const array_1d<double, 3> relaxed_acceleration = (1 - mAlphaBossak) * CurrentAcceleration
  //                                                             + mAlphaBossak * OldAcceleration;
  //         (itNode)->SetValue(RELAXED_ACCELERATION, relaxed_acceleration);
  //     }

  //     //for (ModelPart::ElementsContainerType::ptr_iterator itElem = rModelPart.Elements().ptr_begin(); itElem != rModelPart.Elements().ptr_end(); ++itElem)

  //     #pragma omp parallel for firstprivate(EquationId,RHS_Contribution,LHS_Contribution)
  //     for(int k = 0; k<static_cast<int>(rModelPart.Elements().size()); k++)
  //     {
  //         auto itElem = rModelPart.Elements().ptr_begin()+k;
  //         int thread_id = OpenMPUtils::ThisThread();

  //         (*itElem)->InitializeNonLinearIteration(CurrentProcessInfo);
  //         //KRATOS_WATCH(LHS_Contribution);
  //         //basic operations for the element considered
  //         (*itElem)->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

  //         //std::cout << rCurrentElement->Id() << " RHS = " << RHS_Contribution << std::endl;
  //         (*itElem)->CalculateMassMatrix(mMass[thread_id], CurrentProcessInfo);
  //         (*itElem)->CalculateLocalVelocityContribution(rLHS_Contribution, RHS_Contribution, CurrentProcessInfo);

  //         (*itElem)->EquationIdVector(EquationId, CurrentProcessInfo);

  //         //adding the dynamic contributions (statics is already included)
  //         AddDynamicsToLHS(LHS_Contribution, mDamp[thread_id], mMass[thread_id], CurrentProcessInfo);
  //         AddDynamicsToRHS((*itElem), RHS_Contribution, mDamp[thread_id], mMass[thread_id], CurrentProcessInfo);

  //         Element::GeometryType& rGeom = (*itElem)->GetGeometry();
  //         unsigned int NumNodes = rGeom.PointsNumber();
  //         unsigned int Dimension = rGeom.WorkingSpaceDimension();

  //         unsigned int index = 0;
  //         for (unsigned int i = 0; i < NumNodes; i++)
  //         {
  //             auto& reaction = rGeom[i].FastGetSolutionStepValue(REACTION);

  //             double& target_value0 = reaction[0];
  //             const double& origin_value0 = RHS_Contribution[index++];
  //             #pragma omp atomic
  //             target_value0 -= origin_value0;

  //             double& target_value1 = reaction[1];
  //             const double& origin_value1 = RHS_Contribution[index++];
  //             #pragma omp atomic
  //             target_value1 -= origin_value1;

  //             if (Dimension == 3)
  //             {
  //               double& target_value2 = reaction[2];
  //               const double& origin_value2 = RHS_Contribution[index++];
  //               #pragma omp atomic
  //               target_value2 -= origin_value2;
  //             }
  //     //        rGeom[i].FastGetSolutionStepValue(REACTION_X,0) -= RHS_Contribution[index++];
  //      //          rGeom[i].FastGetSolutionStepValue(REACTION_Y,0) -= RHS_Contribution[index++];
  //     //        if (Dimension == 3) rGeom[i].FastGetSolutionStepValue(REACTION_Z,0) -= RHS_Contribution[index++];
  //             index++; // skip pressure dof
  //         }
  //     }

  //     rModelPart.GetCommunicator().AssembleCurrentData(REACTION);

  //     // Base scheme calls FinalizeSolutionStep method of elements and conditions
  //     Scheme<TSparseSpace, TDenseSpace>::FinalizeSolutionStep(rModelPart, A, Dx, b);
  // }

  //************************************************************************************************
  //************************************************************************************************


  ///@}
  ///@name Inquiry
  ///@{

  ///@}
  ///@name Friends
  ///@{

  ///@}
 protected:
  ///@name Protected static Member Variables
  ///@{

  ///@}
  ///@name Protected member Variables
  ///@{

  ///@}
  ///@name Protected Operators
  ///@{

  ///@}
  ///@name Protected Operations
  ///@{

  /** On periodic boundaries, the nodal area and the values to project need to take into account contributions from elements on
   * both sides of the boundary. This is done using the conditions and the non-historical nodal data containers as follows:\n
   * 1- The partition that owns the PeriodicCondition adds the values on both nodes to their non-historical containers.\n
   * 2- The non-historical containers are added across processes, communicating the right value from the condition owner to all partitions.\n
   * 3- The value on all periodic nodes is replaced by the one received in step 2.
   */
  // void PeriodicConditionProjectionCorrection(ModelPart& rModelPart)
  // {
  //     const int num_nodes = rModelPart.NumberOfNodes();
  //     const int num_conditions = rModelPart.NumberOfConditions();

  //     #pragma omp parallel for
  //     for (int i = 0; i < num_nodes; i++) {
  //         auto it_node = rModelPart.NodesBegin() + i;

  //         it_node->SetValue(NODAL_AREA,0.0);
  //         it_node->SetValue(ADVPROJ,ZeroVector(3));
  //         it_node->SetValue(DIVPROJ,0.0);
  //     }

  //     #pragma omp parallel for
  //     for (int i = 0; i < num_conditions; i++) {
  //         auto it_cond = rModelPart.ConditionsBegin() + i;

  //         if(it_cond->Is(PERIODIC)) {
  //             this->AssemblePeriodicContributionToProjections(it_cond->GetGeometry());
  //         }
  //     }

  //     rModelPart.GetCommunicator().AssembleNonHistoricalData(NODAL_AREA);
  //     rModelPart.GetCommunicator().AssembleNonHistoricalData(ADVPROJ);
  //     rModelPart.GetCommunicator().AssembleNonHistoricalData(DIVPROJ);

  //     #pragma omp parallel for
  //     for (int i = 0; i < num_nodes; i++) {
  //         auto it_node = rModelPart.NodesBegin() + i;
  //         this->CorrectContributionsOnPeriodicNode(*it_node);
  //     }
  // }

  // void AssemblePeriodicContributionToProjections(Geometry< Node<3> >& rGeometry)
  // {
  //     unsigned int nodes_in_cond = rGeometry.PointsNumber();

  //     double nodal_area = 0.0;
  //     array_1d<double,3> momentum_projection = ZeroVector(3);
  //     double mass_projection = 0.0;
  //     for ( unsigned int i = 0; i < nodes_in_cond; i++ )
  //     {
  //         auto& r_node = rGeometry[i];
  //         nodal_area += r_node.FastGetSolutionStepValue(NODAL_AREA);
  //         noalias(momentum_projection) += r_node.FastGetSolutionStepValue(ADVPROJ);
  //         mass_projection += r_node.FastGetSolutionStepValue(DIVPROJ);
  //     }

  //     for ( unsigned int i = 0; i < nodes_in_cond; i++ )
  //     {
  //         auto& r_node = rGeometry[i];
  //         /* Note that this loop is expected to be threadsafe in normal conditions,
  //         * since each node should belong to a single periodic link. However, I am
  //         * setting the locks for openmp in case that we try more complicated things
  //         * in the future (like having different periodic conditions for different
  //         * coordinate directions).
  //         */
  //         r_node.SetLock();
  //         r_node.GetValue(NODAL_AREA) = nodal_area;
  //         noalias(r_node.GetValue(ADVPROJ)) = momentum_projection;
  //         r_node.GetValue(DIVPROJ) = mass_projection;
  //         r_node.UnSetLock();
  //     }
  // }

  // void CorrectContributionsOnPeriodicNode(Node<3>& rNode)
  // {
  //     if (rNode.GetValue(NODAL_AREA) != 0.0) // Only periodic nodes will have a non-historical NODAL_AREA set.
  //     {
  //         rNode.FastGetSolutionStepValue(NODAL_AREA) = rNode.GetValue(NODAL_AREA);
  //         noalias(rNode.FastGetSolutionStepValue(ADVPROJ)) = rNode.GetValue(ADVPROJ);
  //         rNode.FastGetSolutionStepValue(DIVPROJ) = rNode.GetValue(DIVPROJ);
  //     }
  // }

  ///@}
  ///@name Protected  Access
  ///@{

  ///@}
  ///@name Protected Inquiry
  ///@{

  ///@}
  ///@name Protected LifeCycle
  ///@{

  ///@}
 private:
  ///@name Static Member Variables
  ///@{

  ///@}
  ///@name Member Variables
  ///@{

  SlipCoordinateTransformation<LocalSystemMatrixType,LocalSystemVectorType,double> mRotationTool;

  ///@}
  ///@name Private Operators
  ///@{

  ///@}
  ///@name Private Operations
  ///@{

  ///@}
  ///@name Private  Access
  ///@{

  ///@}
  ///@name Serialization
  ///@{

  ///@}
  ///@name Private Inquiry
  ///@{

  ///@}
  ///@name Un accessible methods
  ///@{

  ///@}

}; // Class AleSolutionScheme
///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_ALE_SOLUTION_SCHEME_H_INCLUDED  defined
