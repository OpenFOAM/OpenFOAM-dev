/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "AtomizationModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::AtomizationModel<CloudType>::AtomizationModel
(
    CloudType& owner
)
:
    CloudSubModelBase<CloudType>(owner)
{}


template<class CloudType>
Foam::AtomizationModel<CloudType>::AtomizationModel
(
    const AtomizationModel<CloudType>& am
)
:
    CloudSubModelBase<CloudType>(am)
{}


template<class CloudType>
Foam::AtomizationModel<CloudType>::AtomizationModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::AtomizationModel<CloudType>::~AtomizationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::AtomizationModel<CloudType>::initLiquidCore() const
{
    notImplemented
    (
        "Foam::scalar "
        "Foam::AtomizationModel<CloudType>::initLiquidCore() const"
    );

    return 0.0;
}


template<class CloudType>
Foam::scalar Foam::AtomizationModel<CloudType>::Taverage
(
    const scalar& Tl,
    const scalar& Tc
) const
{
    return (2.0*Tl + Tc)/3.0;
}


template<class CloudType>
bool Foam::AtomizationModel<CloudType>::calcChi() const
{
    notImplemented("bool Foam::AtomizationModel<CloudType>::calcChi()");

    return false;
}


template<class CloudType>
void Foam::AtomizationModel<CloudType>::update
(
    const scalar dt,
    scalar& d,
    scalar& liquidCore,
    scalar& tc,
    const scalar rho,
    const scalar mu,
    const scalar sigma,
    const scalar volFlowRate,
    const scalar rhoAv,
    const scalar Urel,
    const vector& pos,
    const vector& injectionPos,
    const scalar pAmbient,
    const scalar chi,
    cachedRandom& rndGen
) const
{
    notImplemented
    (
        "void Foam::AtomizationModel<CloudType>::update"
        "("
            "const scalar, "
            "scalar&, "
            "scalar&, "
            "scalar&, "
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "const vector&, "
            "const vector&, "
            "const scalar, "
            "const scalar, "
            "cachedRandom&"
        ") const"
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AtomizationModelNew.C"

// ************************************************************************* //

