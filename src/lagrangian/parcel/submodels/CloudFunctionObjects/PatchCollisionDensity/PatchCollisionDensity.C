/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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
#include "nonConformalFvPatch.H"
#include "fvPatchFieldMapper.H"
#include "setSizeFvPatchFieldMapper.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::PatchCollisionDensity<CloudType>::write()
{
    const scalarField z(this->owner().mesh().nCells(), 0);

    volScalarField
    (
        IOobject
        (
            this->owner().name() + ":numberCollisionDensity",
            this->owner().mesh().time().name(),
            this->owner().mesh()
        ),
        this->owner().mesh(),
        dimless/dimArea,
        z,
        numberCollisionDensity_
    )
   .write();

    volScalarField
    (
        IOobject
        (
            this->owner().name() + ":numberCollisionDensityRate",
            this->owner().mesh().time().name(),
            this->owner().mesh()
        ),
        this->owner().mesh(),
        dimless/dimArea/dimTime,
        z,
        (numberCollisionDensity_ - numberCollisionDensity0_)
       /(this->owner().mesh().time().value() - time0_)
    )
   .write();

    volScalarField
    (
        IOobject
        (
            this->owner().name() + ":massCollisionDensity",
            this->owner().mesh().time().name(),
            this->owner().mesh()
        ),
        this->owner().mesh(),
        dimMass/dimArea,
        z,
        massCollisionDensity_
    )
   .write();

    volScalarField
    (
        IOobject
        (
            this->owner().name() + ":massCollisionDensityRate",
            this->owner().mesh().time().name(),
            this->owner().mesh()
        ),
        this->owner().mesh(),
        dimMass/dimArea/dimTime,
        z,
        (massCollisionDensity_ - massCollisionDensity0_)
       /(this->owner().mesh().time().value() - time0_)
    )
   .write();

    numberCollisionDensity0_ == numberCollisionDensity_;
    massCollisionDensity0_ == massCollisionDensity_;
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
    numberCollisionDensity_
    (
        this->owner().mesh().boundary(),
        volScalarField::Internal::null(),
        calculatedFvPatchField<scalar>::typeName
    ),
    numberCollisionDensity0_
    (
        this->owner().mesh().boundary(),
        volScalarField::Internal::null(),
        calculatedFvPatchField<scalar>::typeName
    ),
    massCollisionDensity_
    (
        this->owner().mesh().boundary(),
        volScalarField::Internal::null(),
        calculatedFvPatchField<scalar>::typeName
    ),
    massCollisionDensity0_
    (
        this->owner().mesh().boundary(),
        volScalarField::Internal::null(),
        calculatedFvPatchField<scalar>::typeName
    ),
    time0_(this->owner().mesh().time().value())
{
    numberCollisionDensity_ == 0;
    numberCollisionDensity0_ == 0;
    massCollisionDensity_ == 0;
    massCollisionDensity0_ == 0;

    typeIOobject<volScalarField> numberIo
    (
        this->owner().name() + ":numberCollisionDensity",
        this->owner().mesh().time().name(),
        this->owner().mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (numberIo.headerOk())
    {
        const volScalarField numberCollisionDensity
        (
            numberIo,
            this->owner().mesh()
        );
        numberCollisionDensity_ == numberCollisionDensity.boundaryField();
        numberCollisionDensity0_ == numberCollisionDensity.boundaryField();
    }

    typeIOobject<volScalarField> massIo
    (
        this->owner().name() + ":massCollisionDensity",
        this->owner().mesh().time().name(),
        this->owner().mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (massIo.headerOk())
    {
        const volScalarField massCollisionDensity
        (
            massIo,
            this->owner().mesh()
        );
        massCollisionDensity_ == massCollisionDensity.boundaryField();
        massCollisionDensity0_ == massCollisionDensity.boundaryField();
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
    numberCollisionDensity_
    (
        volScalarField::Internal::null(),
        ppm.numberCollisionDensity_
    ),
    numberCollisionDensity0_
    (
        volScalarField::Internal::null(),
        ppm.numberCollisionDensity0_
    ),
    massCollisionDensity_
    (
        volScalarField::Internal::null(),
        ppm.massCollisionDensity_
    ),
    massCollisionDensity0_
    (
        volScalarField::Internal::null(),
        ppm.massCollisionDensity0_
    ),
    time0_(ppm.time0_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchCollisionDensity<CloudType>::~PatchCollisionDensity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PatchCollisionDensity<CloudType>::preEvolve()
{
    const fvMesh& mesh = this->owner().mesh();

    if (!mesh.conformal())
    {
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& fvp = mesh.boundary()[patchi];

            if (isA<nonConformalFvPatch>(fvp))
            {
                const setSizeFvPatchFieldMapper mapper(fvp.size());

                numberCollisionDensity_[patchi].map
                (
                    numberCollisionDensity_[patchi],
                    mapper
                );
                numberCollisionDensity0_[patchi].map
                (
                    numberCollisionDensity0_[patchi],
                    mapper
                );
                massCollisionDensity_[patchi].map
                (
                    massCollisionDensity_[patchi],
                    mapper
                );
                massCollisionDensity0_[patchi].map
                (
                    massCollisionDensity0_[patchi],
                    mapper
                );

                numberCollisionDensity_[patchi] == 0;
                numberCollisionDensity0_[patchi] == 0;
                massCollisionDensity_[patchi] == 0;
                massCollisionDensity0_[patchi] == 0;
            }
        }
    }
}


template<class CloudType>
void Foam::PatchCollisionDensity<CloudType>::postPatch
(
    const parcelType& p,
    const polyPatch& pp
)
{
    if (pp.coupled()) return;

    const label patchi = pp.index();
    const label patchFacei = p.face() - pp.start();

    vector nw, Up;
    this->owner().patchData(p, pp, nw, Up);

    const scalar speed = (p.U() - Up) & nw;

    if (speed > minSpeed_)
    {
        const scalar magSf =
            this->owner().mesh().magSf().boundaryField()[patchi][patchFacei];

        numberCollisionDensity_[patchi][patchFacei] +=
            p.nParticle()/magSf;
        massCollisionDensity_[patchi][patchFacei] +=
            p.nParticle()*p.mass()/magSf;
    }
}


// ************************************************************************* //
