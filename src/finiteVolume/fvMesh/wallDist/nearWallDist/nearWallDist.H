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
    Foam::nearWallDist

Description
    Distance calculation for cells with face on a wall.
    Searches pointNeighbours to find closest.

SourceFiles
    nearWallDist.C

\*---------------------------------------------------------------------------*/

#ifndef nearWallDist_H
#define nearWallDist_H

#include "DemandDrivenMeshObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class nearWallDist Declaration
\*---------------------------------------------------------------------------*/

class nearWallDist
:
    public DemandDrivenMeshObject
    <
        fvMesh,
        DeletableMeshObject,
        nearWallDist
    >
{
    // Private Data

        //- Near-well distance-to-wall boundary field
        volScalarField::Boundary y_;


protected:

    friend class DemandDrivenMeshObject
    <
        fvMesh,
        DeletableMeshObject,
        nearWallDist
    >;

    // Protected Constructors

        //- Construct from mesh
        explicit nearWallDist(const fvMesh& mesh);


public:

    // Declare name of the class and its debug switch
    ClassName("nearWallDist");


    // Constructors

        //- Disallow default bitwise copy construction
        nearWallDist(const nearWallDist&) = delete;


    //- Destructor
    virtual ~nearWallDist();


    // Member Functions

        //- Access to the near-wall distances
        const volScalarField::Boundary& y() const
        {
            return y_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const nearWallDist&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
