/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::multiSolidBodyMotionSolver

Description
    Solid-body motion of the mesh specified by a run-time selectable
    motion function.

SourceFiles
    multiSolidBodyMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef multiSolidBodyMotionSolver_H
#define multiSolidBodyMotionSolver_H

#include "points0MotionSolver.H"
#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class multiSolidBodyMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class multiSolidBodyMotionSolver
:
    public points0MotionSolver
{
    // Private Data

        //- The motion control function
        PtrList<solidBodyMotionFunction> SBMFs_;

        //- Specified cellZones
        labelList zoneIndices_;

        //- Points to move per cellZone
        labelListList pointIndices_;


public:

    //- Runtime type information
    TypeName("multiSolidBodyMotionSolver");


    // Constructors

        //- Construct from mesh and dictionary
        multiSolidBodyMotionSolver
        (
            const word& name,
            const polyMesh&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        multiSolidBodyMotionSolver(const multiSolidBodyMotionSolver&) = delete;


    //- Destructor
    ~multiSolidBodyMotionSolver();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- This is a solid body motion
        virtual bool solidBody() const
        {
            return true;
        }

        //- Solve for motion
        virtual void solve()
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiSolidBodyMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
