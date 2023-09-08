/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "VolumeFraction.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::VolumeFraction<CloudType>::write()
{
    alpha_.write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VolumeFraction<CloudType>::VolumeFraction
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    alpha_
    (
        IOobject
        (
            this->owner().name() + ":alpha",
            this->owner().mesh().time().name(),
            this->owner().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->owner().mesh(),
        dimensionedScalar(dimless, 0)
    )
{}


template<class CloudType>
Foam::VolumeFraction<CloudType>::VolumeFraction
(
    const VolumeFraction<CloudType>& vf
)
:
    CloudFunctionObject<CloudType>(vf),
    alpha_(this->owner().name() + ":alpha", vf.alpha_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VolumeFraction<CloudType>::~VolumeFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::VolumeFraction<CloudType>::postEvolve()
{
    alpha_ = this->owner().alpha();

    CloudFunctionObject<CloudType>::postEvolve();
}


// ************************************************************************* //
