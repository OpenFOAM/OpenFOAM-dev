/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2026 OpenFOAM Foundation
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
    Foam::pointMeshMovers::linearDisplacement

Description
    Mesh motion solver simple linear expansion and contraction of a mesh
    region defined by a motion axis and the extents of the motion.

    Example:
    \verbatim
        mover           pointMeshMover;

        libs            ("libfvMotionSolvers.so");

        pointMeshMover  linearDisplacement;

        axis            (0 1 0);
        xFixed          0.8;
        xMoving         0;
        displacement    table
        (
            (0 0)
            (4 0.7)
        );
    \endverbatim

    This mesh is compressed between \c xFixed and \c xMoving in the direction
    \c axis between time 0 and 4 with a maximum displacement of 0.7.
    The mesh beyond \c xFixed is fixed and beyond \c xMoving moves with maximum
    displacement.

SourceFiles
    linearDisplacement_pointMeshMover.C

\*---------------------------------------------------------------------------*/

#ifndef linearDisplacement_pointMeshMover_H
#define linearDisplacement_pointMeshMover_H

#include "displacementPoints0.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{

/*---------------------------------------------------------------------------*\
                   Class linearDisplacement Declaration
\*---------------------------------------------------------------------------*/

class linearDisplacement
:
    public displacementPoints0
{
    // Private Data

        //- Axis of the linear motion
        vector axis_;

        //- The location of the end of motion
        scalar xFixed_;

        //- The location moving with the specified displacement
        scalar xMoving_;

        //- The displacement as a function of time
        autoPtr<Function1<scalar>> displacement_;


public:

    //- Runtime type information
    TypeName("linearDisplacement");


    // Constructors

        //- Construct from polyMesh and dictionary
        linearDisplacement(const polyMesh&, const dictionary&);

        //- Disallow default bitwise copy construction
        linearDisplacement
        (
            const linearDisplacement&
        ) = delete;


    //- Destructor
    ~linearDisplacement();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> newPoints();

        //- Update from another mesh using the given map (not implemented)
        virtual void mapMesh(const polyMeshMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const linearDisplacement&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pointMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
