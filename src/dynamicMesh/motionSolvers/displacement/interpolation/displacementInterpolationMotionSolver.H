/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::displacementInterpolationMotionSolver

Description
    Mesh motion solver for a polyMesh.

    Scales in between motion prescribed on faceZones. Works out per point
    the distance between the bounding face zones (in all three directions)
    at the start and then every time step
    - moves the faceZones based on tables
    - interpolates the displacement of all points based on the
      faceZone motion.

    Tables are in the \a constant/tables directory.

SourceFiles
    displacementInterpolationMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementInterpolationMotionSolver_H
#define displacementInterpolationMotionSolver_H

#include "points0MotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class displacementInterpolationMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class displacementInterpolationMotionSolver
:
    public points0MotionSolver
{
    // Private Data

        // Face zone information (note: could pack these to only contain
        // used zones)

            //- Interpolation table. From faceZone to times.
            List<scalarField> times_;

            //- Interpolation table. From faceZone to displacements.
            List<vectorField> displacements_;

        // Range information.

            //- Per direction, per range the index of the lower
            //  faceZone
            FixedList<labelList, 3> rangeToZone_;

            //- Per direction, per range the points that are in it
            FixedList<labelListList, 3> rangeToPoints_;

            //- Per direction, per range the weight of the points relative
            //  to this and the next range.
            FixedList<List<scalarField>, 3> rangeToWeights_;


public:

    //- Runtime type information
    TypeName("displacementInterpolation");


    // Constructors

        //- Construct from polyMesh and dictionary
        displacementInterpolationMotionSolver
        (
            const polyMesh&,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        displacementInterpolationMotionSolver
        (
            const displacementInterpolationMotionSolver&
        );


    //- Destructor
    ~displacementInterpolationMotionSolver();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve()
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const displacementInterpolationMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
