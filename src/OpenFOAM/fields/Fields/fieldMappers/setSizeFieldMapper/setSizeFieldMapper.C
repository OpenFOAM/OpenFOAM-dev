/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "setSizeFieldMapper.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setSizeFieldMapper::setSizeFieldMapper(const label size)
:
    size_(size)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::setSizeFieldMapper::operator()
(
    Field<scalar>& f,
    const Field<scalar>&
) const
{
    setSize(f);
}


void Foam::setSizeFieldMapper::operator()
(
    Field<vector>& f,
    const Field<vector>&
) const
{
    setSize(f);
}


void Foam::setSizeFieldMapper::operator()
(
    Field<sphericalTensor>& f,
    const Field<sphericalTensor>&
) const
{
    setSize(f);
}


void Foam::setSizeFieldMapper::operator()
(
    Field<symmTensor>& f,
    const Field<symmTensor>&
) const
{
    setSize(f);
}


void Foam::setSizeFieldMapper::operator()
(
    Field<tensor>& f,
    const Field<tensor>&
) const
{
    setSize(f);
}


Foam::tmp<Foam::Field<Foam::scalar>> Foam::setSizeFieldMapper::operator()
(
    const Field<scalar>&
) const
{
    return setSize<scalar>();
}


Foam::tmp<Foam::Field<Foam::vector>> Foam::setSizeFieldMapper::operator()
(
    const Field<vector>&
) const
{
    return setSize<vector>();
}


Foam::tmp<Foam::Field<Foam::sphericalTensor>>
Foam::setSizeFieldMapper::operator()
(
    const Field<sphericalTensor>&
) const
{
    return setSize<sphericalTensor>();
}


Foam::tmp<Foam::Field<Foam::symmTensor>> Foam::setSizeFieldMapper::operator()
(
    const Field<symmTensor>&
) const
{
    return setSize<symmTensor>();
}


Foam::tmp<Foam::Field<Foam::tensor>> Foam::setSizeFieldMapper::operator()
(
    const Field<tensor>&
) const
{
    return setSize<tensor>();
}


// ************************************************************************* //
