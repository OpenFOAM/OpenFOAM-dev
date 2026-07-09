/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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
    Foam::solvers::NAME

Description
    Solver module for ...

    Reference:
    \verbatim
        Greenshields, C. J., & Weller, H. G. (2022).
        Notes on Computational Fluid Dynamics: General Principles.
        CFD Direct Ltd.: Reading, UK.
    \endverbatim

SourceFiles
    NAME.C

See also
    Foam::solvers::PARENT

\*---------------------------------------------------------------------------*/

#ifndef NAME_H
#define NAME_H

#include "PARENT.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class NAME Declaration
\*---------------------------------------------------------------------------*/

class NAME
:
    public PARENT
{

protected:

    // Protected data

        //- Velocity field
        volVectorField U_;


    // Protected Member Functions

        //- Example function
        //void functionName();


public:

    // Public Data

        //- Reference to the velocity field
        const volVectorField& U;


    //- Runtime type information
    TypeName("NAME");


    // Constructors

        //- Construct from region mesh
        NAME(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        NAME(const NAME&) = delete;


    //- Destructor
    virtual ~NAME();


    // Member Functions

        //- Called at the start of the time-step, before the PIMPLE loop
        //virtual void preSolve();

        //- Called at the start of the PIMPLE loop to move the mesh
        //virtual void moveMesh();

        //- Corrections that follow mesh motion
        //virtual void motionCorrector();

        //- Called at the start of the PIMPLE loop
        //virtual void prePredictor();

        //- Predict the momentum transport
        //virtual void momentumTransportPredictor();

        //- Predict thermophysical transport
        //virtual void thermophysicalTransportPredictor();

        //- Construct/solve the momentum equation
        //virtual void momentumPredictor();

        //- Construct/solve energy equation, update T and thermophysical model
        //virtual void thermophysicalPredictor();

        //- Construct/solve the pressure equation in the PISO loop
        //virtual void pressureCorrector();

        //- Correct the momentum transport
        //virtual void momentumTransportCorrector();

        //- Correct the thermophysical transport
        //virtual void thermophysicalTransportCorrector();

        //- Called after the PIMPLE loop at the end of the time-step
        //virtual void postSolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const NAME&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
