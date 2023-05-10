/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "Flux.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

template<class CloudType>
const Foam::dimensionSet Foam::NumberFlux<CloudType>::dimensions =
    dimless/dimTime;

template<class CloudType>
const Foam::dimensionSet Foam::VolumeFlux<CloudType>::dimensions =
    dimVolume/dimTime;

template<class CloudType>
const Foam::dimensionSet Foam::MassFlux<CloudType>::dimensions =
    dimMass/dimTime;


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType, class Derived>
void Foam::Flux<CloudType, Derived>::write()
{
    if (write_)
    {
        phi_.write();
    }
}


template<class CloudType, class Derived>
void Foam::Flux<CloudType, Derived>::accumulate
(
    const parcelType& p,
    const bool isPre
)
{
    const polyMesh& mesh = this->owner().mesh();
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    const bool isInternal = p.onInternalFace(mesh);

    const label facei = p.face();
    const label bFacei = isInternal ? -1 : facei - mesh.nInternalFaces();
    const label patchi = isInternal ? -1 : bMesh.patchID()[bFacei];
    const label patchFacei = isInternal ? -1 : bMesh.patchFaceID()[bFacei];

    const bool own = mesh.faceOwner()[facei] == p.cell();

    const scalar sign = own == isPre ? +1 : -1;

    scalar& phif =
        isInternal
      ? phi_[facei]
      : phi_.boundaryFieldRef()[patchi][patchFacei];

    phif += sign*Derived::dPhiDeltaT(p)/mesh.time().deltaTValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType, class Derived>
Foam::Flux<CloudType, Derived>::Flux
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, Derived::typeName),
    write_(dict.lookupOrDefault<bool>("write", true)),
    phi_
    (
        IOobject
        (
            this->owner().name() + ":" + Derived::typeName,
            this->owner().mesh().time().name(),
            this->owner().mesh()
        ),
        this->owner().mesh(),
        dimensionedScalar(Derived::dimensions, 0)
    )
{}


template<class CloudType, class Derived>
Foam::Flux<CloudType, Derived>::Flux
(
    const Flux<CloudType, Derived>& ppm
)
:
    CloudFunctionObject<CloudType>(ppm),
    write_(ppm.write_),
    phi_(ppm.phi_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType, class Derived>
Foam::Flux<CloudType, Derived>::~Flux()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType, class Derived>
void Foam::Flux<CloudType, Derived>::preEvolve()
{
    phi_ = dimensionedScalar(Derived::dimensions, 0);
}


template<class CloudType, class Derived>
void Foam::Flux<CloudType, Derived>::preFace(const parcelType& p)
{
    accumulate(p, true);
}


template<class CloudType, class Derived>
void Foam::Flux<CloudType, Derived>::postFace(const parcelType& p)
{
    if (p.onBoundaryFace(this->owner().mesh()))
    {
        accumulate(p, false);
    }
}


// ************************************************************************* //
