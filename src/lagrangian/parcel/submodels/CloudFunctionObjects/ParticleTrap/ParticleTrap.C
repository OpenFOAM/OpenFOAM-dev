/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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

#include "ParticleTrap.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleTrap<CloudType>::ParticleTrap
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    alphaName_
    (
        this->coeffDict().template lookupOrDefault<word>("alpha", "alpha")
    ),
    alphaPtr_(nullptr),
    gradAlphaPtr_(nullptr),
    threshold_(this->coeffDict().template lookup<scalar>("threshold"))
{}


template<class CloudType>
Foam::ParticleTrap<CloudType>::ParticleTrap
(
    const ParticleTrap<CloudType>& pt
)
:
    CloudFunctionObject<CloudType>(pt),
    alphaName_(pt.alphaName_),
    alphaPtr_(pt.alphaPtr_),
    gradAlphaPtr_(nullptr),
    threshold_(pt.threshold_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleTrap<CloudType>::~ParticleTrap()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleTrap<CloudType>::preEvolve()
{
    if (alphaPtr_ == nullptr)
    {
        const fvMesh& mesh = this->owner().mesh();
        const volScalarField& alpha =
            mesh.lookupObject<volScalarField>(alphaName_);

        alphaPtr_ = &alpha;
    }

    if (gradAlphaPtr_.valid())
    {
        gradAlphaPtr_() == fvc::grad(*alphaPtr_);
    }
    else
    {
        gradAlphaPtr_.reset(new volVectorField(fvc::grad(*alphaPtr_)));
    }
}


template<class CloudType>
void Foam::ParticleTrap<CloudType>::postEvolve()
{
    gradAlphaPtr_.clear();
}


template<class CloudType>
void Foam::ParticleTrap<CloudType>::postMove
(
    parcelType& p,
    const scalar,
    const point&,
    bool&
)
{
    if (alphaPtr_->primitiveField()[p.cell()] < threshold_)
    {
        const vector& gradAlpha = gradAlphaPtr_()[p.cell()];
        vector nHat = gradAlpha/mag(gradAlpha);
        scalar nHatU = nHat & p.U();

        if (nHatU < 0)
        {
            p.U() -= 2*nHat*nHatU;
        }
    }
}


// ************************************************************************* //
