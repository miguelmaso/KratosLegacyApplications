/*
==============================================================================
KratosIncompressibleFluidApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//
//   Project Name:        Kratos
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-10-13 08:17:41 $
//   Revision:            $Revision: 1.20 $
//
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/calculate_forces_utilities.h"
#include "custom_utilities/level_set_utilities.h"
#include "custom_utilities/level_set_utilities_implicitextrapolation.h"
#include "custom_utilities/edge_data.h"
#include "custom_utilities/fluid_solver.h"
// #include "custom_utilities/levelset_fluid_solver.h"
#include "custom_utilities/Turbolence_Smagorinsky.h"
#include "custom_utilities/pure_convection_edgebased.h"
#include "custom_utilities/coupled_eulerian_ulf_utilities.h"
#include "custom_utilities/elembased_distance_utilities.h"
#include "custom_utilities/elembased_extrapolation_utilities.h"
#include "custom_utilities/elembased_BC_utilities.h"
#include "custom_utilities/assign_no_slip_condition.h"
#include "custom_utilities/mark_for_refinement.h"
#include "custom_utilities/parallel_extrapolation_utilities.h"
#include "custom_utilities/wave_generator.h"
#include "custom_utilities/estimate_dt_utilities.h"
#include "custom_utilities/lagrangian_particle_utilities.h"
#include "custom_utilities/embedded_utilities.h"
#include "custom_utilities/fluid_thermal_solver_utilities.h"
#include "custom_utilities/move_particle_utility.h"
#include "custom_utilities/particle_utilities.h"
#include "custom_utilities/combustion_utilities.h"
#include "custom_utilities/lagrangian_pfem2_utilities.h"
//#include "custom_utilities/edgebased_levelset.h"



#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"


namespace Kratos
{

namespace Python
{

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<CoupledEulerianUlfUtils>(m, "CoupledEulerianUlfUtils")
    .def(py::init<>())
    .def("SavePseudoLagPart", &CoupledEulerianUlfUtils::SavePseudoLagPart)
    .def("ApplyProjDirichlet", &CoupledEulerianUlfUtils::ApplyProjDirichlet)
    .def("FindInterface", &CoupledEulerianUlfUtils::FindInterface)
    .def("FindIntersectionOfEdges", &CoupledEulerianUlfUtils::FindIntersectionOfEdges)
    .def("DisableSubdomain", &CoupledEulerianUlfUtils::DisableSubdomain)
    ;

    py::class_<EmbeddedUtils > (m, "EmbeddedUtils")
    .def(py::init<>())
    // .def("SaveInterfaceElemsModelPart", &EmbeddedUtils::SaveInterfaceElemsModelPart)
    //.def "CreateIntersectionConditions", &EmbeddedUtils::CreateIntersectionConditions)
    .def("DisableSubdomain3D", &EmbeddedUtils::DisableSubdomain)
    .def("CreateIntersConditions", &EmbeddedUtils::CreateIntersConditions)
    .def("ApplyProjDirichlet", &EmbeddedUtils::ApplyProjDirichlet)
    ;

    py::class_<CalculateForcesUtils > (m, "CalculateForcesUtils")
    .def(py::init<>())
    .def("CalculateForces3D", &CalculateForcesUtils::CalculateForces3D)
    .def("CalculateForces2D", &CalculateForcesUtils::CalculateForces2D)
    .def("CalculatePressureForces2D", &CalculateForcesUtils::CalculatePressureForces2D)
    ;

    py::class_<RefinementUtilities > (m, "RefinementUtilities")
    .def(py::init<>())
    .def("MarkForRefinement", &RefinementUtilities::MarkForRefinement)
    .def("UpdateErrorRatio", &RefinementUtilities::UpdateErrorRatio)
    .def("RelativeSubscaleErrorEstimate", &RefinementUtilities::SubscaleErrorEstimate)
    ;

    py::class_<AssignNoSlipCondition > (m, "AssignNoSlipCondition")
    .def(py::init<>())
    .def("AssignNoSlipCondition2D", &AssignNoSlipCondition::AssignNoSlipCondition2D)
    ;

    py::class_<SmagorinskyTurbulentModel> (m, "SmagorinskyTurbulentModel")
    .def(py::init<>())
    .def("CalculateTurbulentViscosity2D", &SmagorinskyTurbulentModel::CalculateTurbulentViscosity < 2 >)
    .def("CalculateTurbulentViscosity3D", &SmagorinskyTurbulentModel::CalculateTurbulentViscosity < 3 >)
    ;

    py::class_<LevelSetUtilities > (m, "LevelSetUtilities")
    .def(py::init<>())
    .def("RegenerateFluidModelPart", &LevelSetUtilities::RegenerateFluidModelPart)
    .def("MarkNodesAsVisited", &LevelSetUtilities::MarkNodesAsVisited)
    .def("SetDistanceToNegative", &LevelSetUtilities::SetDistanceToNegative)
    .def("ExtrapolateVelocities", &LevelSetUtilities::ExtrapolateVelocities)
    .def("ExtrapolateVelocitiesByLayer", &LevelSetUtilities::ExtrapolateVelocitiesByLayer)
    .def("GenerateModelPart", &LevelSetUtilities::GenerateModelPart)
    .def("ImplicitExtrapolation_PreProcess", &LevelSetUtilities::ImplicitExtrapolation_PreProcess)
    .def("ImplicitExtrapolation_PostProcess", &LevelSetUtilities::ImplicitExtrapolation_PostProcess)
    .def("PrepareForInternalFluidDistanceComputation", &LevelSetUtilities::PrepareForInternalFluidDistanceComputation)
    .def("FluidDistanceComputation_FromBoundary", &LevelSetUtilities::FluidDistanceComputation_FromBoundary)
    .def("ApplyMinimumExtrapolationPressureFix", &LevelSetUtilities::ApplyMinimumExtrapolationPressureFix)
    ;

    py::class_<LevelSetUtilitiesImplicitExtrapolation > (m, "LevelSetUtilitiesImplicitExtrapolation")
    .def(py::init<>())
    .def("RegenerateFluidModelPart", &LevelSetUtilitiesImplicitExtrapolation::RegenerateFluidModelPart)
    .def("MarkNodesAsVisited", &LevelSetUtilitiesImplicitExtrapolation::MarkNodesAsVisited)
    .def("SetDistanceToNegative", &LevelSetUtilitiesImplicitExtrapolation::SetDistanceToNegative)
    .def("ExtrapolateVelocities", &LevelSetUtilitiesImplicitExtrapolation::ExtrapolateVelocities)
    .def("PrepareForInternalFluidDistanceComputation", &LevelSetUtilitiesImplicitExtrapolation::PrepareForInternalFluidDistanceComputation)
    .def("FluidDistanceComputation_FromBoundary", &LevelSetUtilitiesImplicitExtrapolation::FluidDistanceComputation_FromBoundary)
    .def("ApplyFluidProperties", &LevelSetUtilitiesImplicitExtrapolation::ApplyFluidProperties)
    ;

    py::class_< ElemBasedDistanceUtilities > (m, "ElemBasedDistanceUtilities")
    .def(py::init<ModelPart& >())
    .def("IdentifyFreeSurface", &ElemBasedDistanceUtilities::IdentifyFreeSurface)
    .def("MarkExternalAndMixedNodes", &ElemBasedDistanceUtilities::MarkExternalAndMixedNodes)
    .def("MarkInternalAndMixedNodes", &ElemBasedDistanceUtilities::MarkInternalAndMixedNodes)
    .def("SaveScalarVariableToOldStep", &ElemBasedDistanceUtilities::SaveScalarVariableToOldStep)
    .def("ChangeSignToDistance", &ElemBasedDistanceUtilities::ChangeSignToDistance)
    .def("MarkNodesByDistance", &ElemBasedDistanceUtilities::MarkNodesByDistance)
    ;

    py::class_< ElemBasedExtrapolationUtilities > (m, "ElemBasedExtrapolationUtilities")
    .def(py::init<ModelPart& >())
    .def("ExtrapolateVelocities", &ElemBasedExtrapolationUtilities::ExtrapolateVelocities)
    ;

    py::class_< ElemBasedBCUtilities > (m, "ElemBasedBCUtilities")
    .def(py::init<ModelPart& >())
    .def("SetDividedElem_2D", &ElemBasedBCUtilities::SetDividedElem_2D)
    .def("SetPressureAndVelocityFixities", &ElemBasedBCUtilities::SetPressureAndVelocityFixities)
    .def("FreePressureAndVelocity", &ElemBasedBCUtilities::FreePressureAndVelocity)
    .def("SetToZeroPressureAndVelocity", &ElemBasedBCUtilities::SetToZeroPressureAndVelocity)
    ;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    // 	class_< ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>,  boost::noncopyable >       ("ElementBasedNavierStokesSolver2D", init<ModelPart&>() )
    //                           .def("ConstructSystemStructure",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::ConstructSystemStructure)
    //                           .def("Clear",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::Clear)
    //                           .def("SolveStep1",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::SolveStep1)
    //                           .def("SolveStep2",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::SolveStep2)
    //                           .def("SolveStep3",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::SolveStep3)
    //                         ;
    py::class_< MatrixContainer < 2, SparseSpaceType>, typename MatrixContainer < 2, SparseSpaceType>::Pointer > (m, "MatrixContainer2D")
    .def(py::init< >())
    .def("ConstructCSRVector", &MatrixContainer < 2, SparseSpaceType >::ConstructCSRVector)
    .def("BuildCSRData", &MatrixContainer < 2, SparseSpaceType >::BuildCSRData)
    .def("Clear", &MatrixContainer < 2, SparseSpaceType >::Clear)
    ;

    py::class_< MatrixContainer < 3, SparseSpaceType>, typename MatrixContainer < 3, SparseSpaceType>::Pointer > (m, "MatrixContainer3D")
    .def(py::init< >())
    .def("ConstructCSRVector", &MatrixContainer < 3, SparseSpaceType >::ConstructCSRVector)
    .def("BuildCSRData", &MatrixContainer < 3, SparseSpaceType >::BuildCSRData)
    .def("Clear", &MatrixContainer < 3, SparseSpaceType >::Clear)
    ;

    py::class_<
        FluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType>,
        typename FluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType>::Pointer >
    (m, "FluidSolver2D")
    .def(py::init < MatrixContainer < 2, SparseSpaceType>&, ModelPart&, const double, const double, const Vector, bool, double, double, double, double, bool >())
    .def("Initialize", &FluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
    .def("ComputeTimeStep", &FluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
    .def("SolveStep1", &FluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep1)
    .def("SolveStep2", &FluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep2)
    .def("SolveStep3", &FluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep3)
    .def("ComputeTimeStep", &FluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
    .def("CalculateNormals", &FluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculateNormals)
    .def("UpdateFixedVelocityValues", &FluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::UpdateFixedVelocityValues)
    .def("ComputePressureStabilization", &FluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputePressureStabilization)
    .def("ViscosityCorrectionStep", &FluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ViscosityCorrectionStep)
    .def("ComputeViscousForces", &FluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeViscousForces)
    .def("ComputeReactions", &FluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeReactions)
    .def("ComputeMinimum_Havg", &FluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeMinimum_Havg)
    ;

    py::class_<
        FluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType>,
        typename FluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType>::Pointer>
    (m, "FluidSolver3D")
    .def(py::init < MatrixContainer < 3, SparseSpaceType>&, ModelPart&, const double, const double, const Vector, bool, double, double, double, double, bool >())
    .def("Initialize", &FluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
    .def("ComputeTimeStep", &FluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
    .def("SolveStep1", &FluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep1)
    .def("SolveStep2", &FluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep2)
    .def("SolveStep3", &FluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep3)
    .def("ComputeTimeStep", &FluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
    .def("CalculateNormals", &FluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculateNormals)
    .def("UpdateFixedVelocityValues", &FluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::UpdateFixedVelocityValues)
    .def("ComputePressureStabilization", &FluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputePressureStabilization)
    .def("ViscosityCorrectionStep", &FluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ViscosityCorrectionStep)
    .def("ComputeViscousForces", &FluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeViscousForces)
    .def("ComputeReactions", &FluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeReactions)
    .def("ComputeMinimum_Havg", &FluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeMinimum_Havg)
    ;


//             class_< LevelSetFluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType>, boost::noncopyable > ("LevelSetFluidSolver2D", init < MatrixContainer < 2, SparseSpaceType>&, ModelPart&, bool, bool >())
//                     .def("Initialize", &LevelSetFluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
//                     .def("ComputeTimeStep", &LevelSetFluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
//                     .def("SolveStep1", &LevelSetFluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep1)
//                     .def("SolveStep2", &LevelSetFluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep2)
//                     .def("SolveStep3", &LevelSetFluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep3)
//                     .def("SolveStep4", &LevelSetFluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep4)
//                     .def("ExtrapolateVelocities", &LevelSetFluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ExtrapolateVelocities)
//                     .def("MarkExternalAndMixedNodes", &LevelSetFluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkExternalAndMixedNodes)
//                     .def("MarkInternalAndMixedNodes", &LevelSetFluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalAndMixedNodes)
//                     .def("SaveScalarVariableToOldStep", &LevelSetFluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SaveScalarVariableToOldStep)
//                     .def("ChangeSignToDistance", &LevelSetFluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ChangeSignToDistance)
//                     .def("MarkNodesByDistance", &LevelSetFluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkNodesByDistance)
//                     .def("CalculateForces", &LevelSetFluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculateForces)
//                     .def("CalculateVariablesDistribution", &LevelSetFluidSolver < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculateVariablesDistribution)
//                     ;
//
//             class_< LevelSetFluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType>, boost::noncopyable > ("LevelSetFluidSolver3D", init < MatrixContainer < 3, SparseSpaceType>&, ModelPart&, bool, bool >())
//                     .def("Initialize", &LevelSetFluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
//                     .def("ComputeTimeStep", &LevelSetFluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
//                     .def("SolveStep1", &LevelSetFluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep1)
//                     .def("SolveStep2", &LevelSetFluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep2)
//                     .def("SolveStep3", &LevelSetFluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep3)
//                     .def("SolveStep4", &LevelSetFluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep4)
//                     .def("ExtrapolateVelocities", &LevelSetFluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ExtrapolateVelocities)
//                     .def("MarkExternalAndMixedNodes", &LevelSetFluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkExternalAndMixedNodes)
//                     .def("MarkInternalAndMixedNodes", &LevelSetFluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalAndMixedNodes)
//                     .def("SaveScalarVariableToOldStep", &LevelSetFluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SaveScalarVariableToOldStep)
//                     .def("ChangeSignToDistance", &LevelSetFluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ChangeSignToDistance)
//                     .def("MarkNodesByDistance", &LevelSetFluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkNodesByDistance)
//                     .def("CalculateForces", &LevelSetFluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculateForces)
//                     .def("CalculateVariablesDistribution", &LevelSetFluidSolver < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculateVariablesDistribution)
//                     ;


    py::class_<
        PureConvectionEdgeBased < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType>,
        typename PureConvectionEdgeBased < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType>::Pointer>
    (m, "PureConvectionEdgeBased2D")
    .def(py::init<MatrixContainer < 2, SparseSpaceType>&, ModelPart& >())
    .def("Initialize", &PureConvectionEdgeBased < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
    .def("Solve", &PureConvectionEdgeBased < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Solve)
    .def("ComputeTimeStep", &PureConvectionEdgeBased < 2, MatrixContainer < 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
    ;


    py::class_<
        PureConvectionEdgeBased < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType>,
        typename PureConvectionEdgeBased < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType>::Pointer>
    (m, "PureConvectionEdgeBased3D")
    .def(py::init<MatrixContainer < 3, SparseSpaceType>&, ModelPart& >())
    .def("Initialize", &PureConvectionEdgeBased < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
    .def("Solve", &PureConvectionEdgeBased < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Solve)
    .def("ComputeTimeStep", &PureConvectionEdgeBased < 3, MatrixContainer < 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
    ;

    py::class_<ParallelExtrapolationUtilities < 3 > > (m, "ParallelExtrapolationUtilities3D")
    .def(py::init<>())
    .def("ExtrapolateVelocity", &ParallelExtrapolationUtilities < 3 > ::ExtrapolateVelocity)
    .def("ExtrapolatePressureProjection", &ParallelExtrapolationUtilities < 3 > ::ExtrapolatePressureProjection)
    .def("AssignFreeSurfacePressure", &ParallelExtrapolationUtilities < 3 > ::AssignFreeSurfacePressure)
	.def("ExtrapolateTemperature", &ParallelExtrapolationUtilities < 3 > ::ExtrapolateTemperature)
    ;

    py::class_< WaveGenerator > (m, "WaveGenerator")
    .def(py::init<>())
    .def("GenerateWaveXYPlane", &WaveGenerator::GenerateWaveXYPlane)
    .def("GenerateVolumeWaveXYPlane", &WaveGenerator::GenerateVolumeWaveXYPlane)
    .def("GenerateComposedVolumeWaveXYPlane", &WaveGenerator::GenerateComposedVolumeWaveXYPlane)
    ;

    py::class_< EstimateDtUtil <2> > (m, "EstimateDt2D")
    .def(py::init<ModelPart&>())
    .def("EstimateDt", &EstimateDtUtil < 2 > ::EstimateDt)
    .def("CalculateLocalCFL", &EstimateDtUtil < 2 > ::CalculateLocalCFL)
    ;

    py::class_< EstimateDtUtil <3> > (m, "EstimateDt3D")
    .def(py::init<ModelPart&>())
    .def("EstimateDt", &EstimateDtUtil < 3 > ::EstimateDt)
    .def("CalculateLocalCFL", &EstimateDtUtil < 3 > ::CalculateLocalCFL)
    ;

    py::class_< LagrangianPFEM2Utilities <3> > (m, "LagrangianPFEM2Utilities3D")
    .def(py::init<>())
    .def("DetectInletAndOutlet", &LagrangianPFEM2Utilities < 3 > ::DetectInletAndOutlet)
    .def("MoveMesh_ForwardEuler",	&LagrangianPFEM2Utilities < 3 > ::MoveMesh_ForwardEuler)
    .def("ActOnInlet" ,            	&LagrangianPFEM2Utilities < 3 > ::ActOnInlet)
    .def("ActOnOutlet",    		&LagrangianPFEM2Utilities < 3 > ::ActOnOutlet)
    .def("MarkOuterNodes", 		&LagrangianPFEM2Utilities < 3 > ::MarkOuterNodes)
    .def("MoveMesh_Streamlines", 	&LagrangianPFEM2Utilities < 3 > ::MoveMesh_Streamlines)
    .def("EraseOuterElements", 		&LagrangianPFEM2Utilities < 3 > ::EraseOuterElements)
     .def("MarkExcessivelyCloseNodes", 		&LagrangianPFEM2Utilities < 3 > ::MarkExcessivelyCloseNodes)
    ;

    // 	class_< ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>,  boost::noncopyable >       ("ElementBasedNavierStokesSolver2D", init<ModelPart&>() )
    //                           .def("ConstructSystemStructure",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::ConstructSystemStructure)
    //                           .def("Clear",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::Clear)
    //                           .def("SolveStep1",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::SolveStep1)
    //                           .def("SolveStep2",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::SolveStep2)
    //                           .def("SolveStep3",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::SolveStep3)
    // 			  .def("CalculateProjection",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::CalculateProjection)
    // 			;

    //	  class_< EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType>,  boost::noncopyable >       ("EdgeBasedLevelSet2D", init< MatrixContainer< 2, SparseSpaceType>&, ModelPart& >() )
    //			  .def("Initialize",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
    //			  .def("ComputeTimeStep",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
    //			  .def("SolveStep1",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep1)
    //			  .def("SolveStep2",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep2)
    //			  .def("SolveStep3",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep3)
    //			  .def("ExtrapolateValues",&EdgeBasedLevelSet< 2,MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ExtrapolateValues)
    //			  .def("MarkExternalAndMixedNodes",&EdgeBasedLevelSet< 2,MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkExternalAndMixedNodes)
    //			  .def("MarkInternalAndMixedNodes",&EdgeBasedLevelSet< 2,MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalAndMixedNodes)
    //			  .def("SaveScalarVariableToOldStep",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SaveScalarVariableToOldStep)
    //			  .def("ChangeSignToDistance",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ChangeSignToDistance)
    //			  .def("MarkNodesByDistance",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkNodesByDistance)
    ////			  .def("CalculateForces",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculateForces)
    //			  ;



    py::class_<LagrangianParticleUtils <2> >(m, "LagrangianUtils2D")
    .def(py::init<>())
    .def("StreamlineMove", &LagrangianParticleUtils < 2 > ::StreamlineMove)
    .def("StreamlineCorrect", &LagrangianParticleUtils < 2 > ::StreamlineCorrect)
    .def("Reseed", &LagrangianParticleUtils < 2 > ::Reseed)
    .def("ReseedEmptyElements", &LagrangianParticleUtils < 2 > ::ReseedEmptyElements)
    .def("VisualizationModelPart", &LagrangianParticleUtils < 2 > ::VisualizationModelPart)
    .def("TransferToEulerianMesh", &LagrangianParticleUtils < 2 > ::TransferToEulerianMesh)
    .def("RestartStep", &LagrangianParticleUtils < 2 > ::RestartStep)
    .def("BackAndForth", &LagrangianParticleUtils < 2 > ::BackAndForth)
    .def("ConvectParticles", &LagrangianParticleUtils < 2 > ::ConvectParticles)
    .def("TransferToEulerianMeshShapeBased", &LagrangianParticleUtils < 2 > ::TransferToEulerianMeshShapeBased)
    ;

    py::class_<ParticleUtils <2> >(m, "ParticleUtils2D")
    .def(py::init<>())
    .def("StreamlineMove", &ParticleUtils < 2 > ::StreamlineMove)
    .def("Back", &ParticleUtils < 2 > ::Back)
    .def("Reseed", &ParticleUtils < 2 > ::Reseed)
    .def("ReseedEmptyElements", &ParticleUtils < 2 > ::ReseedEmptyElements)
    .def("VisualizationModelPart", &ParticleUtils < 2 > ::VisualizationModelPart)
    .def("TransferToEulerianMesh", &ParticleUtils < 2 > ::TransferToEulerianMesh)
//                    .def("RestartStep", &ParticleUtils < 2 > ::RestartStep)
    .def("TransferToEulerianMeshShapeBased", &ParticleUtils < 2 > ::TransferToEulerianMeshShapeBased)
    .def("EstimateTime", &ParticleUtils < 2 > ::EstimateTime)
    .def("aa", &ParticleUtils < 2 > ::aa)
    .def("Density", &ParticleUtils < 2 > ::Density)
    .def("Density1", &ParticleUtils < 2 > ::Density1)
    .def("Back1", &ParticleUtils < 2 > ::Back1)
    .def("StreamlineMove2", &ParticleUtils < 2 > ::StreamlineMove2)

//		    .def("StreamlineMove2", &ParticleUtils < 2 > ::StreamlineMove2)
    ;


    py::class_<CombustionUtilities >(m, "CombustionUtilities")
    .def(py::init<>())
    .def("Mixture_Fraction",&CombustionUtilities::Mixture_Fraction)
    .def("Enthalpy",&CombustionUtilities::Enthalpy)
    .def("Temperature",&CombustionUtilities::Temperature)
    ;


    py::class_<FluidThermalSolverUtilities>(m, "FluidThermalSolverUtilities")
    .def(py::init<ModelPart&, ModelPart&>())
    .def("ProjectFromThermalToFluid",&FluidThermalSolverUtilities::ProjectFromThermalToFluid)
    .def("ProjectFromFluidToThermal",&FluidThermalSolverUtilities::ProjectFromFluidToThermal)
    .def("ApplyTables",&FluidThermalSolverUtilities::ApplyTables)
    ;


    py::class_< MoveParticleUtility<2> > (m, "MoveParticleUtility2D")
    .def(py::init<ModelPart& , ModelPart& >())
    .def("MountBin", &MoveParticleUtility<2>::MountBin)
    .def("MoveParticles", &MoveParticleUtility<2>::MoveParticles)
    ;

}

} // namespace Python.

} // Namespace Kratos

