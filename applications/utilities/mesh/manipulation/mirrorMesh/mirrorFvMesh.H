/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::mirrorFvMesh

Description

SourceFiles
    mirrorFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef mirrorFvMesh_H
#define mirrorFvMesh_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class mirrorFvMesh Declaration
\*---------------------------------------------------------------------------*/

class mirrorFvMesh
:
    public fvMesh
{
    // Private data

        //- Mirror dictionary
        IOdictionary mirrorMeshDict_;

        //- Mirrored mesh
        fvMesh* mirrorMeshPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        mirrorFvMesh(const mirrorFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const mirrorFvMesh&);


public:


    // Constructors

        //- Construct from IOobject
        mirrorFvMesh(const IOobject& io);


    //- Destructor
    ~mirrorFvMesh();


    // Member Functions

        //- Return reference to mirror mesh
        const fvMesh& mirrorMesh() const
        {
            return *mirrorMeshPtr_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
