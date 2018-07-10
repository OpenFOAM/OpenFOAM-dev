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
    Foam::scalarRanges

Description
    A List of scalarRange.

SourceFiles
    scalarRanges.C

\*---------------------------------------------------------------------------*/

#ifndef scalarRanges_H
#define scalarRanges_H

#include "scalarRange.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class scalarRanges Declaration
\*---------------------------------------------------------------------------*/

class scalarRanges
:
    public List<scalarRange>
{
public:

    // Constructors

        //- Construct Null
        scalarRanges();

        //- Construct from Istream.
        //  The list items are comma-delimited.
        scalarRanges(Istream&);

    // Member Functions

        //- Return true if the given value is within the ranges
        bool selected(const scalar) const;

        //- Return the set of selected entries in the given list
        //  that are within the ranges
        List<bool> selected(const List<scalar>&) const;

        //- Select a list of values that are within the ranges
        List<scalar> select(const List<scalar>&) const;

        //- Select a list of values that are within the ranges
        void inplaceSelect(List<scalar>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
