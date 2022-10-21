/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::solvers::fluid

Description
    Solver module for steady or transient turbulent flow of compressible fluids
    with heat-transfer for HVAC and similar applications, with optional
    mesh motion and change.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, Lagrangian particles,
    radiation, surface film etc. and constraining or limiting the solution.

    Reference:
    \verbatim
        Greenshields, C. J., & Weller, H. G. (2022).
        Notes on Computational Fluid Dynamics: General Principles.
        CFD Direct Ltd.: Reading, UK.
    \endverbatim

SourceFiles
    fluid.C

\*---------------------------------------------------------------------------*/

#ifndef fluid_H
#define fluid_H

#include "isothermalFluid.H"
#include "fluidThermoThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class fluid Declaration
\*---------------------------------------------------------------------------*/

class fluid
:
    public isothermalFluid
{

protected:

    // Thermophysical transport

        autoPtr<fluidThermoThermophysicalTransportModel>
            thermophysicalTransport;


public:

    //- Runtime type information
    TypeName("fluid");


    // Constructors

        //- Construct from region mesh
        fluid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        fluid(const fluid&) = delete;


    //- Destructor
    virtual ~fluid();


    // Member Functions

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Correct the thermophysical transport modelling
        virtual void thermophysicalTransportCorrector();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fluid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
