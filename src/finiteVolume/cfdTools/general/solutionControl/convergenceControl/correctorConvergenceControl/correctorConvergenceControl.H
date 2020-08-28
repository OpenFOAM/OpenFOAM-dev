/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 OpenFOAM Foundation
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
    Foam::correctorConvergenceControl

Description
    Corrector convergence control class. Provides methods to check the
    convergence of an inner iteration loop (e.g., pimple) against both absolute
    and relative residual tolerances.

SourceFiles
    correctorConvergenceControl.C
    correctorConvergenceControlTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef correctorConvergenceControl_H
#define correctorConvergenceControl_H

#include "fvMesh.H"
#include "solutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class correctorConvergenceControl Declaration
\*---------------------------------------------------------------------------*/

class correctorConvergenceControl
{
public:

    //- Residual correction data structure
    struct corrResidualData
    {
        wordRe name;
        scalar absTol;
        scalar relTol;
    };


private:

    // Private Member Functions

        //- Disallow default bitwise copy construction
        correctorConvergenceControl
        (
            const correctorConvergenceControl&
        ) = delete;

        //- Disallow default bitwise assignment
        void operator=(const correctorConvergenceControl&) = delete;


protected:

    // Protected data

        //- Reference to the solution control
        const solutionControl& control_;

        //- Name of the corrector loop
        const word loopName_;


    // Protected Member Functions

        //- Get the number of solves that have happened for this variable in
        //  this time-step
        static void getNSolves
        (
            const fvMesh& mesh,
            const word& fieldName,
            label& n
        );


    // Protected Static Template Functions

        //- Get the number of solves that have happened for this variable in
        //  this time-step
        template<class Type>
        static void getNTypeSolves
        (
            const fvMesh& mesh,
            const word& fieldName,
            label& n
        );


public:

    // Static Data Members

        //- Run-time type information
        TypeName("correctorConvergenceControl");


    // Constructors

        //- Construct from a solution control and the loop name
        correctorConvergenceControl
        (
            const solutionControl& control,
            const word& loopName
        );


    //- Destructor
    virtual ~correctorConvergenceControl();


    // Member Functions

        // Evolution

            //- Return true if corrector residual controls are present
            virtual bool hasCorrResidualControls() const = 0;

            //- Return true if all correction convergence checks are satisfied
            virtual bool corrCriteriaSatisfied() const = 0;

            //- Reset the solve index in the correction residual control data
            virtual void resetCorrSolveIndex() = 0;

            //- Update the solve index in the correction residual control data
            virtual void updateCorrSolveIndex() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "correctorConvergenceControlTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
