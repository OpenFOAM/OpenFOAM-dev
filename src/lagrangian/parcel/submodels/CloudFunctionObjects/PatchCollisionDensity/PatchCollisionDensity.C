/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 OpenFOAM Foundation
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

#include "PatchCollisionDensity.H"
#include "Pstream.H"
#include "stringListOps.H"
#include "ListOps.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::PatchCollisionDensity<CloudType>::write()
{
    const scalarField z(this->owner().mesh().nCells(), 0);

    volScalarField
    (
        IOobject
        (
            this->owner().name() + ":collisionDensity",
            this->owner().mesh().time().timeName(),
            this->owner().mesh()
        ),
        this->owner().mesh(),
        dimless/dimArea,
        z,
        collisionDensity_
    )
   .write();

    volScalarField
    (
        IOobject
        (
            this->owner().name() + ":collisionDensityRate",
            this->owner().mesh().time().timeName(),
            this->owner().mesh()
        ),
        this->owner().mesh(),
        dimless/dimArea/dimTime,
        z,
        (collisionDensity_ - collisionDensity0_)
       /(this->owner().mesh().time().value() - time0_)
    )
   .write();

    collisionDensity0_ == collisionDensity_;
    time0_ = this->owner().mesh().time().value();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchCollisionDensity<CloudType>::PatchCollisionDensity
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    minSpeed_(dict.lookupOrDefault<scalar>("minSpeed", -1)),
    collisionDensity_
    (
        this->owner().mesh().boundary(),
        volScalarField::Internal::null(),
        calculatedFvPatchField<scalar>::typeName
    ),
    collisionDensity0_
    (
        this->owner().mesh().boundary(),
        volScalarField::Internal::null(),
        calculatedFvPatchField<scalar>::typeName
    ),
    time0_(this->owner().mesh().time().value())
{
    collisionDensity_ == 0;
    collisionDensity0_ == 0;

    IOobject io
    (
        this->owner().name() + ":collisionDensity",
        this->owner().mesh().time().timeName(),
        this->owner().mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<volScalarField>())
    {
        const volScalarField collisionDensity(io, this->owner().mesh());
        collisionDensity_ == collisionDensity.boundaryField();
        collisionDensity0_ == collisionDensity.boundaryField();
    }
}


template<class CloudType>
Foam::PatchCollisionDensity<CloudType>::PatchCollisionDensity
(
    const PatchCollisionDensity<CloudType>& ppm
)
:
    CloudFunctionObject<CloudType>(ppm),
    minSpeed_(ppm.minSpeed_),
    collisionDensity_
    (
        volScalarField::Internal::null(),
        ppm.collisionDensity_
    ),
    collisionDensity0_
    (
        volScalarField::Internal::null(),
        ppm.collisionDensity0_
    ),
    time0_(ppm.time0_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchCollisionDensity<CloudType>::~PatchCollisionDensity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PatchCollisionDensity<CloudType>::postPatch
(
    const parcelType& p,
    const polyPatch& pp,
    bool&
)
{
    const label patchi = pp.index();
    const label patchFacei = p.face() - pp.start();

    vector nw, Up;
    this->owner().patchData(p, pp, nw, Up);

    const scalar speed = (p.U() - Up) & nw;
    if (speed > minSpeed_)
    {
        collisionDensity_[patchi][patchFacei] +=
            1/this->owner().mesh().magSf().boundaryField()[patchi][patchFacei];
    }
}


// ************************************************************************* //
