// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_NODAL_CONCENTRATED_ELEMENT_H_INCLUDED )
#define  KRATOS_NODAL_CONCENTRATED_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"

namespace Kratos
{
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

/**
 * @class NodalConcentratedElement
 * @ingroup StructuralMechanicsApplication
 * @brief Concentrated nodal for 3D and 2D points
 * @details The element can consider both the displacement and rotational stiffness, and both the mass and the inertia
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) NodalConcentratedElement
    : public Element
{
public:

    ///@name Type Definitions
    ///@{

    /// Definition of the node type
    typedef Node<3> NodeType;

    /// Definition of the geometry
    typedef Geometry<NodeType> GeometryType;

    /// Definition of the base type
    typedef Element BaseType;

    /// Definition of the index type
    typedef std::size_t IndexType;

    /// Definition of the size type
    typedef std::size_t SizeType;

    /// Counted pointer of NodalConcentratedElement
    KRATOS_CLASS_POINTER_DEFINITION( NodalConcentratedElement);

     /**
     * @brief Flags related to the element computation
     */
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_DISPLACEMENT_STIFFNESS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_NODAL_MASS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_ROTATIONAL_STIFFNESS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_NODAL_INERTIA );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_DAMPING_RATIO );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_ROTATIONAL_DAMPING_RATIO );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RAYLEIGH_DAMPING );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_ACTIVE_NODE_FLAG );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_ONLY_COMPRESSION );

    ///@}
    ///@name  Enum's
    ///@{

    /**
     * @brief If the elemnt works only in compression/tension or both
     */
    enum class CompressionTension {COMPRESSION_AND_TENSION = 0, COMPRESSION = 1, TENSION = 2};

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructors
     * @param NewId The Id of the new created element
     * @param pGeometry Pointer to the node geometry that defines the element
     * @param UseRayleighDamping If use Raleigh damping or damping ratio
     * @param ComputeActiveNodeFlag If the element activates or deactivates in function of the node that define it
     * @param ComputeCompressionTension If the elements works only in compression or tension or both
     */
    NodalConcentratedElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        const bool UseRayleighDamping = false,
        const bool ComputeActiveNodeFlag = true,
        const CompressionTension ComputeCompressionTension = CompressionTension::COMPRESSION_AND_TENSION
        );

    /**
     * @brief Default constructors with properties
     * @param NewId The Id of the new created element
     * @param pGeometry Pointer to the node geometry that defines the element
     * @param pProperties The properties of the element
     * @param UseRayleighDamping If use Raleigh damping or damping ratio
     * @param ComputeActiveNodeFlag If the element activates or deactivates in function of the node that define it
     * @param ComputeCompressionTension If the elements works only in compression or tension or both
     */
    NodalConcentratedElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        const bool UseRayleighDamping = false,
        const bool ComputeActiveNodeFlag = true,
        const CompressionTension ComputeCompressionTension = CompressionTension::COMPRESSION_AND_TENSION
        );

    ///Copy constructor
    NodalConcentratedElement(NodalConcentratedElement const& rOther);

    /// Destructor.
    ~NodalConcentratedElement() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    NodalConcentratedElement& operator=(NodalConcentratedElement const& rOther);

    ///@}
    ///@name Operations
    ///@{
    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    //************* GETTING METHODS

    /**
     * Sets on rElementalDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(
        DofsVectorType& rElementalDofList, 
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     */
    void EquationIdVector(
        EquationIdVectorType& rResult, 
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Sets on rValues the nodal displacements
     */
    void GetValuesVector(Vector& rValues, int Step = 0) override;

    /**
     * Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) override;

    /**
     * Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override;

    //************* STARTING - ENDING  METHODS

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    void Initialize() override;

    /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Called at the end of eahc solution step
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;


    //************* COMPUTING  METHODS

    /**
     * This is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix, 
        VectorType& rRightHandSideVector, 
        ProcessInfo& rCurrentProcessInfo
        ) override;


    /**
     * This calculates just the RHS
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */

    void CalculateRightHandSide( 
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * This calculates just the LHS
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */

    void CalculateLeftHandSide( 
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * this is called during the assembling process in order
      * to calculate the elemental mass matrix
      * @param rMassMatrix: the elemental mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix, 
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * this is called during the assembling process in order
      * to calculate the elemental damping matrix
      * @param rDampingMatrix: the elemental damping matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix, 
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{
    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
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

    Flags mELementalFlags; /// Elemental flags

    ///@name Protected Operators
    ///@{
    NodalConcentratedElement() : Element()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method computes the actual size of the system of equations
     * @return This method returns the size of the system of equations
     */
    std::size_t ComputeSizeOfSystem();

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
    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class NodalConcentratedElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_NODAL_CONCENTRATED_ELEMENT_H_INCLUDED  defined 
