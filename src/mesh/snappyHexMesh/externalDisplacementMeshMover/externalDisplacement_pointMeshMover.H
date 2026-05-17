/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2026 OpenFOAM Foundation
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
    Foam::pointMeshMovers::externalDisplacement

Description
    Mesh motion solver for an fvMesh.  Based on solving the cell-centre
    Laplacian for the motion displacement.

SourceFiles
    externalDisplacement_pointMeshMover.C

\*---------------------------------------------------------------------------*/

#ifndef externalDisplacement_pointMeshMover_H
#define externalDisplacement_pointMeshMover_H

#include "displacement_pointMeshMover.H"
#include "externalDisplacementMeshMover.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{

/*---------------------------------------------------------------------------*\
                  Class externalDisplacement Declaration
\*---------------------------------------------------------------------------*/

class externalDisplacement
:
    public pointMeshMovers::displacement
{
    // Private Data

        dictionary dict_;

        mutable autoPtr<externalDisplacementMeshMover> meshMoverPtr_;


public:

    //- Runtime type information
    TypeName("externalDisplacement");


    // Constructors

        //- Construct from polyMesh and dictionary
        externalDisplacement(const polyMesh&, const dictionary&);

        //- Disallow default bitwise copy construction
        externalDisplacement(const externalDisplacement&) = delete;


    //- Destructor
    ~externalDisplacement();


    // Member Functions

        externalDisplacementMeshMover& meshMover() const;

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> newPoints();

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //- Update topology
        virtual void topoChange(const polyTopoChangeMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const externalDisplacement&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pointMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
