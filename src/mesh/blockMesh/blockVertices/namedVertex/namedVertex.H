/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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
    Foam::blockVertices::namedVertex

Description
    Gives name to a vertex.

SourceFiles
    namedVertex.C

\*---------------------------------------------------------------------------*/

#ifndef blockVertices_namedVertex_H
#define blockVertices_namedVertex_H

#include "blockVertex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockVertices
{

/*---------------------------------------------------------------------------*\
                         Class namedVertex Declaration
\*---------------------------------------------------------------------------*/

class namedVertex
:
    public blockVertex
{
protected:

    // Protected member data

        //- The dictionary variable name for the vertex number
        word name_;

        //- The vertex location
        autoPtr<blockVertex> vertexPtr_;


public:

    //- Runtime type information
    TypeName("name");


    // Constructors

        //- Construct from Istream setting pointsList
        namedVertex
        (
            const dictionary&,
            const label index,
            const searchableSurfaces& geometry,
            Istream&
        );


    //- Destructor
    virtual ~namedVertex()
    {}


    // Member Functions

        virtual operator point() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockVertices
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
