/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "SinglePhaseMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::SinglePhaseMixture<CloudType>::constructIds()
{
    if (this->phaseProps().size() == 0)
    {
        FatalErrorInFunction
            << "Phase list is empty" << exit(FatalError);
    }
    else if (this->phaseProps().size() > 1)
    {
        FatalErrorInFunction
            << "Only one phase permitted" << exit(FatalError);
    }

    switch (this->phaseProps()[0].phase())
    {
        case phaseProperties::GAS:
        {
            idGas_ = 0;
            break;
        }
        case phaseProperties::LIQUID:
        {
            idLiquid_ = 0;
            break;
        }
        case phaseProperties::SOLID:
        {
            idSolid_ = 0;
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown phase enumeration" << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SinglePhaseMixture<CloudType>::SinglePhaseMixture
(
    const dictionary& dict,
    CloudType& owner
)
:
    CompositionModel<CloudType>(dict, owner, typeName),
    idGas_(-1),
    idLiquid_(-1),
    idSolid_(-1)
{
    constructIds();
}


template<class CloudType>
Foam::SinglePhaseMixture<CloudType>::SinglePhaseMixture
(
    const SinglePhaseMixture<CloudType>& cm
)
:
    CompositionModel<CloudType>(cm),
    idGas_(cm.idGas_),
    idLiquid_(cm.idLiquid_),
    idSolid_(cm.idSolid_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SinglePhaseMixture<CloudType>::~SinglePhaseMixture()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const Foam::scalarField&
Foam::SinglePhaseMixture<CloudType>::YMixture0() const
{
    return this->phaseProps()[0].Y();
}


template<class CloudType>
Foam::label Foam::SinglePhaseMixture<CloudType>::idGas() const
{
    return idGas_;
}


template<class CloudType>
Foam::label Foam::SinglePhaseMixture<CloudType>::idLiquid() const
{
    return idLiquid_;
}


template<class CloudType>
Foam::label Foam::SinglePhaseMixture<CloudType>::idSolid() const
{
    return idSolid_;
}


// ************************************************************************* //
