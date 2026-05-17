/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2026 OpenFOAM Foundation
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
    Foam::fvMotionSolver

Description
    Base class for fvMesh motion solvers.

SourceFiles
    fvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef fvMotionSolver_H
#define fvMotionSolver_H

#include "fvMeshMover.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class fvMotionSolver
:
    public fvMeshMover
{
protected:

    // Protected Member Functions

        //- Create the corresponding patch types for cellMotion from those
        //  of the given pointMotion
        template<class Type>
        wordList cellMotionBoundaryTypes
        (
            const typename PointField<Type>::
            Boundary& pmUbf
        ) const;


public:

    //- Runtime type information
    ClassName("fvMotionSolver");


    // Constructors

        //- Construct from fvMesh
        fvMotionSolver(fvMesh&);

        //- Construct from polyMesh
        fvMotionSolver(const polyMesh&);


    // Member Functions

        //- Is this motion solid body? Delegate to the motion solver.
        virtual bool solidBodyMotion() const
        {
            return false;
        }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> newPoints() = 0;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        // //- Update corresponding to the given map
        // virtual void topoChange(const polyTopoChangeMap&);

        // //- Update from another mesh using the given map
        // virtual void mapMesh(const polyMeshMap&);

        // //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMotionSolverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
