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
//   Last modified by:    $Author: kazem $
//   Date:                $Date: 2009-01-15 18:46:07 $
//   Revision:            $Revision: 1.15 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "incompressible_fluid_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_io_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_edgebased_levelset_solver_to_python.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;



PYBIND11_MODULE(KratosIncompressibleFluidApplication, m)
{

    py::class_<KratosIncompressibleFluidApplication,
        KratosIncompressibleFluidApplication::Pointer,
        KratosApplication>(m, "KratosIncompressibleFluidApplication")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomIOToPython(m);
    AddCustomProcessesToPython(m);
    AddCustomEdgeBasedLevelSetToPython(m);

    //registering variables in python
//		KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(FRACT_VEL)
    //	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(PRESS_PROJ)
//		KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(CONV_PROJ)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ACTIVATE_TAU2)


    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MACH_NUMBER)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, REACTION_WATER_PRESSURE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BDF_COEFFICIENTS);
//		KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_MASS)
//		KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUX_INDEX)
//		KRATOS_REGISTER_IN_PYTHON_VARIABLE( EXTERNAL_PRESSURE)
// 		KRATOS_REGISTER_IN_PYTHON_VARIABLE( DIAMETER)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PERMEABILITY_INV)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WATER_PRESSURE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AIR_PRESSURE)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AIR_SOUND_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WATER_SOUND_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOUND_VELOCITY)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MIN_DT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MAX_DT)
    //for disabling parts of Eulerian domain
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DISABLE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INLET_VELOCITY)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INLET_PRESSURE)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRESSURE_DT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WATER_PRESSURE_DT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AIR_PRESSURE_DT)

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, RHS_VECTOR)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, AUX_VECTOR)
    //KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(AUX_VEL)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, CONVECTION_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, AUX_VEL1)

//     KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_DIVIDED)
  }

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
