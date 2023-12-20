/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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
    Foam::solidBodyMotionSolver

Description
    Solid-body motion of the mesh specified by a run-time selectable
    motion function.

SourceFiles
    solidBodyMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef solidBodyMotionSolver_H
#define solidBodyMotionSolver_H

#include "points0MotionSolver.H"
#include "solidBodyMotionFunction.H"
#include "polyCellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class solidBodyMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class solidBodyMotionSolver
:
    public points0MotionSolver
{
    // Private Data

        //- The motion control function
        autoPtr<solidBodyMotionFunction> SBMFPtr_;

        //- The name of the cell zone to restrict motion to
        polyCellSet set_;

        //- Points to move when cell zone is supplied
        labelList setPointIndices_;

        //- Cache the current transform for topoChange
        mutable septernion transform_;


    // Private Member Functions

        //- Re-calculate the set point ID-s from the set, following mesh change
        void updateSetPointIndices();


public:

    //- Runtime type information
    TypeName("solidBody");


    // Constructors

        //- Construct from mesh and dictionary
        solidBodyMotionSolver
        (
            const word& name,
            const polyMesh&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        solidBodyMotionSolver(const solidBodyMotionSolver&) = delete;


    //- Destructor
    ~solidBodyMotionSolver();


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

        //- Update local data for topology changes
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solidBodyMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
