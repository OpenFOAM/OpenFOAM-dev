/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "VoidFraction.H"

// * * * * * * * * * * * * * Protectd Member Functions * * * * * * * * * * * //

template<class CloudType>
void Foam::VoidFraction<CloudType>::write()
{
    if (thetaPtr_.valid())
    {
        thetaPtr_->write();
    }
    else
    {
        FatalErrorIn("void Foam::VoidFraction<CloudType>::write()")
            << "thetaPtr not valid" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VoidFraction<CloudType>::VoidFraction
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    thetaPtr_(NULL)
{}


template<class CloudType>
Foam::VoidFraction<CloudType>::VoidFraction
(
    const VoidFraction<CloudType>& vf
)
:
    CloudFunctionObject<CloudType>(vf),
    thetaPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VoidFraction<CloudType>::~VoidFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::VoidFraction<CloudType>::preEvolve()
{
    if (thetaPtr_.valid())
    {
        thetaPtr_->internalField() = 0.0;
    }
    else
    {
        const fvMesh& mesh = this->owner().mesh();

        thetaPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "Theta",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }
}


template<class CloudType>
void Foam::VoidFraction<CloudType>::postEvolve()
{
    volScalarField& theta = thetaPtr_();

    const fvMesh& mesh = this->owner().mesh();

    theta.internalField() /= mesh.time().deltaTValue()*mesh.V();

    CloudFunctionObject<CloudType>::postEvolve();
}


template<class CloudType>
void Foam::VoidFraction<CloudType>::postMove
(
    parcelType& p,
    const label cellI,
    const scalar dt,
    const point&,
    bool&
)
{
    volScalarField& theta = thetaPtr_();

    theta[cellI] += dt*p.nParticle()*p.volume();
}


// ************************************************************************* //
