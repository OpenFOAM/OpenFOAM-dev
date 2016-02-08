/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "Function1.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::TableBase<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os  << nl << indent << table_ << token::END_STATEMENT << nl;
    writeEntries(os);
}


template<class Type>
void Foam::Function1Types::TableBase<Type>::writeEntries(Ostream& os) const
{
    if (boundsHandling_ != CLAMP)
    {
        os.writeKeyword("outOfBounds") << boundsHandlingToWord(boundsHandling_)
            << token::END_STATEMENT << nl;
    }
    if (interpolationScheme_ != "linear")
    {
        os.writeKeyword("interpolationScheme") << interpolationScheme_
            << token::END_STATEMENT << nl;
    }
}


// ************************************************************************* //
