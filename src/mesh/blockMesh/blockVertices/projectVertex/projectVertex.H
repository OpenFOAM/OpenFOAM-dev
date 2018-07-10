/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::blockVertices::projectVertex

Description
    Projects the vertex onto the selected surfacees of the
    geometry provided as a searchableSurfaces object.

SourceFiles
    projectVertex.C

\*---------------------------------------------------------------------------*/

#ifndef blockVertices_projectVertex_H
#define blockVertices_projectVertex_H

#include "pointVertex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockVertices
{

/*---------------------------------------------------------------------------*\
                        Class projectVertex Declaration
\*---------------------------------------------------------------------------*/

class projectVertex
:
    public pointVertex
{
    // Private data

        const searchableSurfaces& geometry_;

        //- The indices of surfaces onto which the points are projected
        labelList surfaces_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        projectVertex(const projectVertex&);

        //- Disallow default bitwise assignment
        void operator=(const projectVertex&);


public:

    //- Runtime type information
    TypeName("project");


    // Constructors

        //- Construct from Istream setting pointsList
        projectVertex
        (
            const dictionary&,
            const label index,
            const searchableSurfaces& geometry,
            Istream&
        );


    //- Destructor
    virtual ~projectVertex()
    {}


    // Member Functions

        //- Project the given points onto the surface
        virtual operator point() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockVertices
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
