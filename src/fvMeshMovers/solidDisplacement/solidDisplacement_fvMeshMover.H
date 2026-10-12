/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::fvMeshMovers::solidDisplacement

Description
    Mesh mover which interpolates the displacement \c volVectorField D

    from the mesh cell and boundary face centres to points using \c
    volPointInterpolation and moves the mesh accordingly.

Usage
    Example:
    \verbatim
    mover
    {
        type            solidDisplacement;

        libs            ("libsolidDisplacement_fvMeshMover.so");
    }
    \endverbatim

SourceFiles
    solidDisplacementFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef solidDisplacement_fvMeshMovers_H
#define solidDisplacement_fvMeshMovers_H

#include "displacement_pointMeshMover.H"
#include "fvMeshMover.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshMovers
{

/*---------------------------------------------------------------------------*\
             Class solidDisplacement Declaration
\*---------------------------------------------------------------------------*/

class solidDisplacement
:
    public fvMeshMover,
    public pointMeshMovers::displacement
{

public:

    //- Runtime type information
    TypeName("solidDisplacement");


    // Constructors

        //- Construct from fvMesh and dictionary
        solidDisplacement(fvMesh&, const dictionary&);

        //- Disallow default bitwise copy construction
        solidDisplacement(const solidDisplacement&) = delete;


    //- Destructor
    ~solidDisplacement();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> newPoints();

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&)
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solidDisplacement&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pointMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
