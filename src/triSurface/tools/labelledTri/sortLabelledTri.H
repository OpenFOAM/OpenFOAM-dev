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
    Foam::sortLabelledTri

Description
    Helper class which when constructed with a triSurface
    sorts the faces according to region number (or rather constructs a
    mapping).

SourceFiles
    sortLabelledTri.C

\*---------------------------------------------------------------------------*/

#ifndef sortLabelledTri_H
#define sortLabelledTri_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class sortLabelledTri;
class triSurface;

/*---------------------------------------------------------------------------*\
                        Class surfAndLabel Declaration
\*---------------------------------------------------------------------------*/

//- Hold surface and label
class surfAndLabel
{
    const triSurface* surfPtr_;

    label index_;

    // Private Classes

        //- Scalar comparison function used for sorting
        class less
        {
        public:

            inline bool operator()
            (
                const surfAndLabel& one,
                const surfAndLabel& two
            ) const;
        };


public:

    friend class sortLabelledTri;

    // Constructors

        //- Construct null
        surfAndLabel()
        :
            surfPtr_(nullptr),
            index_(-1)
        {}

        //- Construct from surface and index
        surfAndLabel(const triSurface& surf, const label index)
        :
            surfPtr_(&surf),
            index_(index)
        {}
};


/*---------------------------------------------------------------------------*\
                       Class sortLabelledTri Declaration
\*---------------------------------------------------------------------------*/

class sortLabelledTri
:
    public List<surfAndLabel>
{

public:

    // Constructors

        //- Construct from surface, sorting the faces according to patch
        sortLabelledTri(const triSurface&);


    // Member Functions

        // Access

            //- Set the labelList to those of sorted point indices
            void indices(labelList&) const;

            //- Return the list of sorted point indices
            labelList indices() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
