/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::dimensionSet& Foam::geometricOneField::dimensions() const
{
    return dimless;
}

inline Foam::scalar Foam::geometricOneField::operator[](const label) const
{
    return scalar(1);
}

inline Foam::oneField Foam::geometricOneField::field() const
{
    return oneField();
}

inline Foam::oneField Foam::geometricOneField::oldTime() const
{
    return oneField();
}

inline Foam::oneFieldField Foam::geometricOneField::boundaryField() const
{
    return oneFieldField();
}


inline const Foam::geometricOneField& Foam::operator*
(
    const geometricOneField& gof,
    const geometricOneField&
)
{
    return gof;
}

inline const Foam::geometricOneField& Foam::operator/
(
    const geometricOneField& gof,
    const geometricOneField&
)
{
    return gof;
}


// ************************************************************************* //
