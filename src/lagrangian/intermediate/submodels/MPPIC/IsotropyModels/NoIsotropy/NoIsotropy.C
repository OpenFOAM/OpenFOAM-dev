/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "NoIsotropy.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IsotropyModels::NoIsotropy<CloudType>::NoIsotropy
(
    const dictionary& dict,
    CloudType& owner
)
:
    IsotropyModel<CloudType>(owner)
{}


template<class CloudType>
Foam::IsotropyModels::NoIsotropy<CloudType>::NoIsotropy
(
    const NoIsotropy<CloudType>& cm
)
:
    IsotropyModel<CloudType>(cm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IsotropyModels::NoIsotropy<CloudType>::~NoIsotropy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::IsotropyModels::NoIsotropy<CloudType>::calculate()
{}


template<class CloudType>
bool Foam::IsotropyModels::NoIsotropy<CloudType>::active() const
{
    return false;
}


// ************************************************************************* //
