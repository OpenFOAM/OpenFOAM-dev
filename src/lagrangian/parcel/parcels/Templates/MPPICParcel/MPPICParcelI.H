/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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

#include "MPPICParcel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::MPPICParcel<ParcelType>::MPPICParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label facei
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti, facei),
    id_(-1, -1)
{}


template<class ParcelType>
inline Foam::MPPICParcel<ParcelType>::MPPICParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    ParcelType(mesh, position, celli),
    id_(-1, -1)
{}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::labelPair& Foam::MPPICParcel<ParcelType>::id() const
{
    return id_;
}


template<class ParcelType>
inline Foam::labelPair& Foam::MPPICParcel<ParcelType>::id()
{
    return id_;
}


// ************************************************************************* //
