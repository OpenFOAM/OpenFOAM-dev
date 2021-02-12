/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2021 OpenFOAM Foundation
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

#include "DenseDragForce.H"
#include "SchillerNaumannDragForce.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DenseDragForce<CloudType>::DenseDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& typeName
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    alphacName_(this->coeffs().lookup("alphac")),
    alphacPtr_(nullptr),
    alphacInterpPtr_(nullptr)
{}


template<class CloudType>
Foam::DenseDragForce<CloudType>::DenseDragForce
(
    const DenseDragForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df),
    alphacName_(df.alphacName_),
    alphacPtr_(nullptr),
    alphacInterpPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DenseDragForce<CloudType>::~DenseDragForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const Foam::interpolation<Foam::scalar>&
Foam::DenseDragForce<CloudType>::alphacInterp() const
{
    if (!alphacInterpPtr_.valid())
    {
        FatalErrorInFunction
            << "Carrier phase volume-fraction interpolation object not set"
            << abort(FatalError);
    }

    return alphacInterpPtr_();
}


template<class CloudType>
void Foam::DenseDragForce<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        if (!this->mesh().template foundObject<volScalarField>(alphacName_))
        {
            alphacPtr_.reset
            (
                new volScalarField(alphacName_, 1 - this->owner().theta())
            );
        }

        const volScalarField& alphac =
            this->mesh().template lookupObject<volScalarField>(alphacName_);

        alphacInterpPtr_.reset
        (
            interpolation<scalar>::New
            (
                this->owner().solution().interpolationSchemes(),
                alphac
            ).ptr()
        );
    }
    else
    {
        alphacInterpPtr_.clear();
        alphacPtr_.clear();
    }
}


// ************************************************************************* //
