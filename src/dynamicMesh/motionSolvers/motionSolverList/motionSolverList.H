/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2021 OpenFOAM Foundation
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
    Foam::motionSolverList

Description
    Motion of the mesh specified as a list of motion solvers.

    The motion solvers are executed in order and the resulting displacements
    accumulated into an overall displacement and the displaced point positions
    returned.

SourceFiles
    motionSolverList.C

\*---------------------------------------------------------------------------*/

#ifndef motionSolverList_H
#define motionSolverList_H

#include "motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class motionSolverList Declaration
\*---------------------------------------------------------------------------*/

class motionSolverList
:
    public motionSolver
{
    // Private Data

        PtrList<motionSolver> motionSolvers_;


public:

    //- Runtime type information
    TypeName("motionSolverList");


    // Constructors

        //- Construct from mesh and dictionary
        motionSolverList(const polyMesh&, const dictionary&);

        //- Disallow default bitwise copy construction
        motionSolverList(const motionSolverList&) = delete;


    //- Destructor
    virtual ~motionSolverList();


    // Member Functions

        //- Provide current points for motion.  Uses current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //-  Update local data for topology changes
        virtual void updateMesh(const mapPolyMesh&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const mapDistributePolyMesh&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const motionSolverList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
