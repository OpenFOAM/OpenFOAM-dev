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

#include "distributedWeightedFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelListList&
Foam::distributedWeightedFvPatchFieldMapper::addressing() const
{
    return addressing_;
}


const Foam::scalarListList&
Foam::distributedWeightedFvPatchFieldMapper::weights() const
{
    return weights_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::distributedWeightedFvPatchFieldMapper::operator()
(
    const Field<scalar>& mapF
) const
{
    return map(mapF);
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::distributedWeightedFvPatchFieldMapper::operator()
(
    const Field<vector>& mapF
) const
{
    return map(mapF);
}


Foam::tmp<Foam::Field<Foam::sphericalTensor>>
Foam::distributedWeightedFvPatchFieldMapper::operator()
(
    const Field<sphericalTensor>& mapF
) const
{
    return map(mapF);
}


Foam::tmp<Foam::Field<Foam::symmTensor>>
Foam::distributedWeightedFvPatchFieldMapper::operator()
(
    const Field<symmTensor>& mapF
) const
{
    return map(mapF);
}


Foam::tmp<Foam::Field<Foam::tensor>>
Foam::distributedWeightedFvPatchFieldMapper::operator()
(
    const Field<tensor>& mapF
) const
{
    return map(mapF);
}


void Foam::distributedWeightedFvPatchFieldMapper::operator()
(
    Field<scalar>& f,
    const Field<scalar>& mapF
) const
{
    return map(f, mapF);
}


void Foam::distributedWeightedFvPatchFieldMapper::operator()
(
    Field<vector>& f,
    const Field<vector>& mapF
) const
{
    return map(f, mapF);
}


void Foam::distributedWeightedFvPatchFieldMapper::operator()
(
    Field<sphericalTensor>& f,
    const Field<sphericalTensor>& mapF
) const
{
    return map(f, mapF);
}


void Foam::distributedWeightedFvPatchFieldMapper::operator()
(
    Field<symmTensor>& f,
    const Field<symmTensor>& mapF
) const
{
    return map(f, mapF);
}


void Foam::distributedWeightedFvPatchFieldMapper::operator()
(
    Field<tensor>& f,
    const Field<tensor>& mapF
) const
{
    return map(f, mapF);
}


// ************************************************************************* //
