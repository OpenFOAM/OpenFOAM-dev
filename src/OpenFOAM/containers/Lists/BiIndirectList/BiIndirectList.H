/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::BiIndirectList

Description
    Indexes into negList (negative index) or posList (zero or positive index).

SourceFiles
    BiIndirectListI.H

\*---------------------------------------------------------------------------*/

#ifndef BiIndirectList_H
#define BiIndirectList_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class BiIndirectList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class BiIndirectList
{
    // Private data

        UList<T>& posList_;
        UList<T>& negList_;
        List<label> addressing_;


public:

    // Constructors

        //- Construct given the complete lists and the addressing array
        inline BiIndirectList
        (
            const UList<T>& posList,
            const UList<T>& negList,
            const labelUList&
        );

        //- Construct given the complete list and by transferring addressing
        inline BiIndirectList
        (
            const UList<T>& posList,
            const UList<T>& negList,
            const Xfer<List<label>>&
        );


    // Member Functions

        // Access

            //- Return the number of elements in the list
            inline label size() const;

            //- Return true if the list is empty (ie, size() is zero).
            inline bool empty() const;

            inline const UList<T>& posList() const;
            inline const UList<T>& negList() const;

            //- Return the list addressing
            inline const List<label>& addressing() const;

            //- Calculate index given whether index is into posList or negList
            inline static label posIndex(const label);
            inline static label negIndex(const label);

        // Edit

            //- Reset addressing
            inline void resetAddressing(const labelUList&);
            inline void resetAddressing(const Xfer<List<label>>&);


        // Member Operators

            //- Return the addressed elements as a List
            inline List<T> operator()() const;

            //- Return non-const access to an element
            inline T& operator[](const label);

            //- Return const access to an element
            inline const T& operator[](const label) const;

            //- Assignment to UList of addressed elements
            inline void operator=(const UList<T>&);

            //- Assignment of all entries to the given value
            inline void operator=(const T&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BiIndirectListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
