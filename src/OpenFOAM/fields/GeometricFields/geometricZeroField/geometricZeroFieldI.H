/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "geometricZeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::dimensionSet& Foam::geometricZeroField::dimensions() const
{
    return dimless;
}


inline Foam::scalar Foam::geometricZeroField::operator[](const label) const
{
    return scalar(0);
}


inline Foam::geometricZeroField Foam::geometricZeroField::oldTime() const
{
    return geometricZeroField();
}


inline Foam::geometricZeroField::Internal
Foam::geometricZeroField::operator()() const
{
    return zeroField();
}


inline Foam::geometricZeroField::Internal
Foam::geometricZeroField::internalField() const
{
    return zeroField();
}


inline Foam::geometricZeroField::Internal
Foam::geometricZeroField::v() const
{
    return Internal();
}


inline typename Foam::geometricZeroField::Internal::FieldType
Foam::geometricZeroField::primitiveField() const
{
    return zeroField();
}


inline Foam::geometricZeroField::Boundary
Foam::geometricZeroField::boundaryField() const
{
    return zeroFieldField();
}


// ************************************************************************* //
