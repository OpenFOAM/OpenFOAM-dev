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
    Foam::topoChangerFvMesh

Description
    Abstract base class for a topology changing fvMesh

SourceFiles
    topoChangerFvMesh.C
    newTopoFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef topoChangerFvMesh_H
#define topoChangerFvMesh_H

#include "dynamicFvMesh.H"
#include "polyTopoChanger.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class topoChangerFvMesh Declaration
\*---------------------------------------------------------------------------*/

class topoChangerFvMesh
:
    public dynamicFvMesh
{
protected:

    // Protected Data

        polyTopoChanger topoChanger_;


public:

    //- Runtime type information
    TypeName("topoChangerFvMesh");


    // Constructors

        //- Construct from objectRegistry, and read/write options
        explicit topoChangerFvMesh(const IOobject& io);

        //- Disallow default bitwise copy construction
        topoChangerFvMesh(const topoChangerFvMesh&) = delete;


    //- Destructor
    virtual ~topoChangerFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const topoChangerFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
