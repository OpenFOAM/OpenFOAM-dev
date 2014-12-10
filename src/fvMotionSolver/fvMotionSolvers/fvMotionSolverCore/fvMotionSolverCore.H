/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::fvMotionSolverCore

Description
    Base class for fvMesh based motionSolvers.

SourceFiles
    fvMotionSolverCore.C

\*---------------------------------------------------------------------------*/

#ifndef fvMotionSolverCore_H
#define fvMotionSolverCore_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class velocityMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class fvMotionSolverCore
{

protected:

    // Protected data

        //- The fvMesh to be moved
        const fvMesh& fvMesh_;


    // Protected Member Functions

        //- Create the corresponding patch types for cellMotion from those
        //  of the given pointMotion
        template<class Type>
        wordList cellMotionBoundaryTypes
        (
            const typename GeometricField<Type, pointPatchField, pointMesh>::
            GeometricBoundaryField& pmUbf
        ) const;

public:

    //- Runtime type information
    ClassName("displacementMotionSolver");

    // Constructors

        //- Construct from polyMesh
        fvMotionSolverCore(const polyMesh&);


    // Member Functions

        //- Return reference to the fvMesh to be moved
        const fvMesh& mesh() const
        {
            return fvMesh_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvMotionSolverCoreTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
