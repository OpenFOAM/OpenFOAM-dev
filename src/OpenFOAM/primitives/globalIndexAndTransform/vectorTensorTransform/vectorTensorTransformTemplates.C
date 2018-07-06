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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::vectorTensorTransform::transform(const Type& x) const
{
    if (hasR_)
    {
        return Foam::transform(R(), x);
    }
    else
    {
        return x;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::vectorTensorTransform::transform
(
    const Field<Type>& fld
) const
{
    if (hasR_)
    {
        return Foam::transform(R(), fld);
    }
    else
    {
        return fld;
    }
}


template<class Type>
Type Foam::vectorTensorTransform::invTransform(const Type& x) const
{
    if (hasR_)
    {
        return Foam::transform(R().T(), x);
    }
    else
    {
        return x;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::vectorTensorTransform::invTransform
(
    const Field<Type>& fld
) const
{
    if (hasR_)
    {
        return Foam::transform(R().T(), fld);
    }
    else
    {
        return fld;
    }
}


// ************************************************************************* //
