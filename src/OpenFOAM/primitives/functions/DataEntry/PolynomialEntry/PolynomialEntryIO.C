/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "PolynomialEntry.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PolynomialEntry<Type>& poly
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const DataEntry<Type>& >(poly)
            << token::SPACE << poly.coeffs_;
    }
    else
    {
        os  << static_cast<const DataEntry<Type>& >(poly);
        os.write
        (
            reinterpret_cast<const char*>(&poly.coeffs_),
            sizeof(poly.coeffs_)
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const PolynomialEntry&)"
    );

    return os;
}


template<class Type>
void Foam::PolynomialEntry<Type>::writeData(Ostream& os) const
{
    DataEntry<Type>::writeData(os);

    os  << nl << indent << coeffs_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
