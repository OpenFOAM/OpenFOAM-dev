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

#include "generalFieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelUList& Foam::generalFieldMapper::directAddressing() const
{
    FatalErrorInFunction
        << "attempt to access null direct addressing"
        << abort(FatalError);

    return labelUList::null();
}


const Foam::labelListList& Foam::generalFieldMapper::addressing() const
{
    FatalErrorInFunction
        << "attempt to access null interpolation addressing"
        << abort(FatalError);

    return labelListList::null();
}


const Foam::scalarListList& Foam::generalFieldMapper::weights() const
{
    FatalErrorInFunction
        << "attempt to access null interpolation weights"
        << abort(FatalError);

    return scalarListList::null();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::generalFieldMapper::operator()
(
    Field<scalar>& f,
    const Field<scalar>& mapF
) const
{
    map(f, mapF);
}


void Foam::generalFieldMapper::operator()
(
    Field<vector>& f,
    const Field<vector>& mapF
) const
{
    map(f, mapF);
}


void Foam::generalFieldMapper::operator()
(
    Field<sphericalTensor>& f,
    const Field<sphericalTensor>& mapF
) const
{
    map(f, mapF);
}


void Foam::generalFieldMapper::operator()
(
    Field<symmTensor>& f,
    const Field<symmTensor>& mapF
) const
{
    map(f, mapF);
}


void Foam::generalFieldMapper::operator()
(
    Field<tensor>& f,
    const Field<tensor>& mapF
) const
{
    map(f, mapF);
}


Foam::tmp<Foam::Field<Foam::scalar>> Foam::generalFieldMapper::operator()
(
    const Field<scalar>& mapF
) const
{
    return map(mapF);
}


Foam::tmp<Foam::Field<Foam::vector>> Foam::generalFieldMapper::operator()
(
    const Field<vector>& mapF
) const
{
    return map(mapF);
}


Foam::tmp<Foam::Field<Foam::sphericalTensor>>
Foam::generalFieldMapper::operator()
(
    const Field<sphericalTensor>& mapF
) const
{
    return map(mapF);
}


Foam::tmp<Foam::Field<Foam::symmTensor>> Foam::generalFieldMapper::operator()
(
    const Field<symmTensor>& mapF
) const
{
    return map(mapF);
}


Foam::tmp<Foam::Field<Foam::tensor>> Foam::generalFieldMapper::operator()
(
    const Field<tensor>& mapF
) const
{
    return map(mapF);
}


// ************************************************************************* //
