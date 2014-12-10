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
    Foam::staticFvMesh

Description
    Foam::staticFvMesh

SourceFiles
    staticFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef staticFvMesh_H
#define staticFvMesh_H

#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class staticFvMesh Declaration
\*---------------------------------------------------------------------------*/

class staticFvMesh
:
    public dynamicFvMesh
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        staticFvMesh(const staticFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const staticFvMesh&);


public:

    //- Runtime type information
    TypeName("staticFvMesh");


    // Constructors

        //- Construct from IOobject
        staticFvMesh(const IOobject& io);


    //- Destructor
    ~staticFvMesh();


    // Member Functions

        //- Dummy update function which does not change the mesh
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
