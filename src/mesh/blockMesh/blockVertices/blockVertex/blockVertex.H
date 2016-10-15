/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::blockVertex

Description
    Define a block vertex.

SourceFiles
    blockVertex.C

\*---------------------------------------------------------------------------*/

#ifndef blockVertex_H
#define blockVertex_H

#include "searchableSurfaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class blockVertex Declaration
\*---------------------------------------------------------------------------*/

class blockVertex
{

public:

    //- Runtime type information
    TypeName("blockVertex");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            blockVertex,
            Istream,
            (
                const searchableSurfaces& geometry,
                Istream& is
            ),
            (geometry, is)
        );


    // Constructors

        //- Construct null
        blockVertex();

        //- Clone function
        virtual autoPtr<blockVertex> clone() const;

        //- New function which constructs and returns pointer to a blockVertex
        static autoPtr<blockVertex> New
        (
            const searchableSurfaces& geometry,
            Istream&
        );

        //- Class used for the read-construction of
        //  PtrLists of blockVertex
        class iNew
        {
            const searchableSurfaces& geometry_;

        public:

            iNew(const searchableSurfaces& geometry)
            :
                geometry_(geometry)
            {}

            autoPtr<blockVertex> operator()(Istream& is) const
            {
                return blockVertex::New(geometry_, is);
            }
        };


    //- Destructor
    virtual ~blockVertex()
    {}


    // Member Functions

        virtual operator point() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
