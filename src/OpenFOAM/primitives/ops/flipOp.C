/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "flipOp.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
Foam::scalar Foam::flipOp::operator()(const scalar& v) const
{
    return -v;
}


template<> Foam::vector Foam::flipOp::operator()(const vector& v) const
{
    return -v;
}


template<>Foam::sphericalTensor Foam::flipOp::operator()
(
    const sphericalTensor& v
) const
{
    return -v;
}


template<> Foam::symmTensor Foam::flipOp::operator()
(
    const symmTensor& v
) const
{
    return -v;
}


template<> Foam::tensor Foam::flipOp::operator()(const tensor& v) const
{
    return -v;
}


template<> Foam::triad Foam::flipOp::operator()
(
    const triad& v
) const
{
    return -v;
}


// ************************************************************************* //
