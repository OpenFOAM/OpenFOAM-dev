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

#include "fieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::mapDistributeBase& Foam::fieldMapper::distributeMap() const
{
    FatalErrorInFunction
        << "attempt to access null distributeMap"
        << abort(FatalError);
    return *(new mapDistributeBase());
}


const Foam::labelUList& Foam::fieldMapper::directAddressing() const
{
    FatalErrorInFunction
        << "attempt to access null direct addressing"
        << abort(FatalError);

    return labelUList::null();
}


const Foam::labelListList& Foam::fieldMapper::addressing() const
{
    FatalErrorInFunction
        << "attempt to access null interpolation addressing"
        << abort(FatalError);

    return labelListList::null();
}


const Foam::scalarListList& Foam::fieldMapper::weights() const
{
    FatalErrorInFunction
        << "attempt to access null interpolation weights"
        << abort(FatalError);

    return scalarListList::null();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalar>> Foam::fieldMapper::operator()
(
    const Field<scalar>& mapF
) const
{
    return map(mapF);
}


Foam::tmp<Foam::Field<Foam::vector>> Foam::fieldMapper::operator()
(
    const Field<vector>& mapF
) const
{
    return map(mapF);
}


Foam::tmp<Foam::Field<Foam::sphericalTensor>> Foam::fieldMapper::operator()
(
    const Field<sphericalTensor>& mapF
) const
{
    return map(mapF);
}


Foam::tmp<Foam::Field<Foam::symmTensor>> Foam::fieldMapper::operator()
(
    const Field<symmTensor>& mapF
) const
{
    return map(mapF);
}


Foam::tmp<Foam::Field<Foam::tensor>> Foam::fieldMapper::operator()
(
    const Field<tensor>& mapF
) const
{
    return map(mapF);
}


void Foam::fieldMapper::operator()
(
    Field<scalar>& f,
    const Field<scalar>& mapF
) const
{
    return map(f, mapF);
}


void Foam::fieldMapper::operator()
(
    Field<vector>& f,
    const Field<vector>& mapF
) const
{
    return map(f, mapF);
}


void Foam::fieldMapper::operator()
(
    Field<sphericalTensor>& f,
    const Field<sphericalTensor>& mapF
) const
{
    return map(f, mapF);
}


void Foam::fieldMapper::operator()
(
    Field<symmTensor>& f,
    const Field<symmTensor>& mapF
) const
{
    return map(f, mapF);
}


void Foam::fieldMapper::operator()
(
    Field<tensor>& f,
    const Field<tensor>& mapF
) const
{
    return map(f, mapF);
}


// ************************************************************************* //
