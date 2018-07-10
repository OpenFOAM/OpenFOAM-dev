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
    Foam::refineCell

Description
    Container with cells to refine. Refinement given as single direction.

SourceFiles
    refineCell.C

\*---------------------------------------------------------------------------*/

#ifndef refineCell_H
#define refineCell_H

#include "label.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class refineCell;

Ostream& operator<<(Ostream&, const refineCell&);


/*---------------------------------------------------------------------------*\
                           Class refineCell Declaration
\*---------------------------------------------------------------------------*/

class refineCell
{
    // Private data

    //- Cell label
    label cellNo_;

    //- Preferred refinement direction (always normalized).
    vector direction_;

public:

    // Constructors

        //- Null
        refineCell();

        //- From components. Vector will be normalized upon construction.
        refineCell(const label, const vector&);

        //- From Istream. Vector will be normalized upon construction.
        refineCell(Istream& is);


    // Member Functions

        label cellNo() const
        {
            return cellNo_;
        }

        const vector& direction() const
        {
            return direction_;
        }


    // Friend Operators

        inline friend bool operator==
        (
            const refineCell& rc1,
            const refineCell& rc2
        )
        {
            return (rc1.cellNo() == rc2.cellNo());
        }

        inline friend bool operator!=
        (
            const refineCell& rc1,
            const refineCell& rc2
        )
        {
            return !(rc1 == rc2);
        }


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const refineCell&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
