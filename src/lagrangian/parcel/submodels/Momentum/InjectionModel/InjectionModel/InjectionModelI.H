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

#include "InjectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::timeStart() const
{
    return SOI_;
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::volumeTotal() const
{
    return volumeTotal_;
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::massTotal() const
{
    return massTotal_;
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::massInjected() const
{
    return massInjected_;
}


template<class CloudType>
Foam::label Foam::InjectionModel<CloudType>::nInjections() const
{
    return nInjections_;
}


template<class CloudType>
Foam::label Foam::InjectionModel<CloudType>::parcelsAddedTotal() const
{
    return parcelsAddedTotal_;
}


// ************************************************************************* //
