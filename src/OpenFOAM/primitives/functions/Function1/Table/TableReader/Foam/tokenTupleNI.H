/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "tokenTupleN.H"
#include "ITstream.H"
#include "SubList.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::label Foam::tokenTupleN::size() const
{
    return offsets_.size() - 1;
}


template<class Type>
inline Type Foam::tokenTupleN::get(const IOstream& s, const label i)
{
    if (i < 0 || i > size() - 1)
    {
        FatalIOErrorInFunction(s)
            << "index " << i << " out of range 0 ... " << size() - 1
            << exit(FatalIOError);
    }

    ITstream its
    (
        s.name(),
        SubList<token>(tokens_, offsets_[i + 1] - offsets_[i], offsets_[i])
    );

    return pTraits<Type>(its);
}


// ************************************************************************* //
