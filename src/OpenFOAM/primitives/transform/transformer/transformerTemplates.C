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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::transformer::transform(const Type& x) const
{
    if (rotates_)
    {
        return Foam::transform(R(), x);
    }
    else
    {
        return x;
    }
}


template<class Type>
void Foam::transformer::transform
(
    Field<Type>& res,
    const Field<Type>& fld
) const
{
    if (rotates_)
    {
        return Foam::transform(res, R(), fld);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::transformer::transform
(
    const Field<Type>& fld
) const
{
    if (rotates_)
    {
        return Foam::transform(R(), fld);
    }
    else
    {
        return fld;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::transformer::transform
(
    const tmp<Field<Type>>& tfld
) const
{
    if (rotates_)
    {
        return Foam::transform(R(), tfld);
    }
    else
    {
        return tfld;
    }
}


template<class Type>
Type Foam::transformer::invTransform(const Type& x) const
{
    if (rotates_)
    {
        return Foam::transform(R().T(), x);
    }
    else
    {
        return x;
    }
}


template<class Type>
void Foam::transformer::invTransform
(
    Field<Type>& res,
    const Field<Type>& fld
) const
{
    if (rotates_)
    {
        Foam::transform(res, R().T(), fld);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::transformer::invTransform
(
    const Field<Type>& fld
) const
{
    if (rotates_)
    {
        return Foam::transform(R().T(), fld);
    }
    else
    {
        return fld;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::transformer::invTransform
(
    const tmp<Field<Type>>& tfld
) const
{
    if (rotates_)
    {
        return Foam::transform(R().T(), tfld);
    }
    else
    {
        return tfld;
    }
}


// ************************************************************************* //
