/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::solvers::functions

Description
    Solver module to execute the \c functionObjects for a specified solver

    The solver specified by either the \c subSolver or if not present the \c
    solver entry in the \c controlDict is instantiated to provide the physical
    fields needed by the \c functionObjects.  The \c functionObjects are then
    instantiated from the specifications are read from the \c functions entry in
    the \c controlDict and executed in a time-loop also controlled by entries in
    \c controlDict and the \c maxDeltaT() returned by the sub-solver.

    The fields and other objects registered by the sub-solver are set to
    NO_WRITE as they are not changed by the execution of the functionObjects and
    should not be written out each write-time.  Fields and other objects created
    and changed by the execution of the functionObjects are written out.

SourceFiles
    functions.C

\*---------------------------------------------------------------------------*/

#ifndef solvers_functions_H
#define solvers_functions_H

#include "movingMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class functions Declaration
\*---------------------------------------------------------------------------*/

class functions
:
    public movingMesh
{

protected:

    //- Pointer to the solid or fluid solver
    autoPtr<solver> solverPtr;


public:

    //- Runtime type information
    TypeName("functions");


    // Constructors

        //- Construct from region mesh
        functions(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        functions(const functions&) = delete;


    //- Destructor
    virtual ~functions();


    // Member Functions

        //- Return the current maximum time-step for stable solution
        virtual scalar maxDeltaT() const;

        //- Called at the beginning of the PIMPLE loop
        virtual void prePredictor();

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Construct and solve the pressure equation in the PISO loop
        virtual void pressureCorrector();

        //- Correct the thermophysical transport modelling
        virtual void postCorrector();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const functions&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
