/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "ParticleTracks.H"
#include "Pstream.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleTracks<CloudType>::write()
{
    if (cloudPtr_.valid())
    {
        cloudPtr_->write();

        if (resetOnWrite_)
        {
            cloudPtr_->clear();
        }
    }
    else
    {
        if (debug)
        {
            InfoInFunction << "cloupPtr invalid" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleTracks<CloudType>::ParticleTracks
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    trackInterval_(this->coeffDict().template lookup<label>("trackInterval")),
    maxSamples_(this->coeffDict().template lookup<label>("maxSamples")),
    resetOnWrite_(this->coeffDict().lookup("resetOnWrite")),
    faceHitCounter_(),
    cloudPtr_(nullptr)
{}


template<class CloudType>
Foam::ParticleTracks<CloudType>::ParticleTracks
(
    const ParticleTracks<CloudType>& ppm
)
:
    CloudFunctionObject<CloudType>(ppm),
    trackInterval_(ppm.trackInterval_),
    maxSamples_(ppm.maxSamples_),
    resetOnWrite_(ppm.resetOnWrite_),
    faceHitCounter_(ppm.faceHitCounter_),
    cloudPtr_(ppm.cloudPtr_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleTracks<CloudType>::~ParticleTracks()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleTracks<CloudType>::preEvolve()
{
    if (!cloudPtr_.valid())
    {
        cloudPtr_.reset
        (
            this->owner().cloneBare(this->owner().name() + "Tracks").ptr()
        );
    }
}


template<class CloudType>
void Foam::ParticleTracks<CloudType>::postFace(const parcelType& p, bool&)
{
    if
    (
        this->owner().solution().output()
     || this->owner().solution().transient()
    )
    {
        if (!cloudPtr_.valid())
        {
            FatalErrorInFunction
             << "Cloud storage not allocated" << abort(FatalError);
        }

        hitTableType::iterator iter =
            faceHitCounter_.find(labelPair(p.origProc(), p.origId()));

        label localI = -1;
        if (iter != faceHitCounter_.end())
        {
            iter()++;
            localI = iter();
        }
        else
        {
            localI = 1;
            faceHitCounter_.insert(labelPair(p.origProc(), p.origId()), localI);
        }

        label nSamples = floor(localI/trackInterval_);
        if ((localI % trackInterval_ == 0) && (nSamples < maxSamples_))
        {
            cloudPtr_->append
            (
                static_cast<parcelType*>(p.clone(this->owner().mesh()).ptr())
            );
        }
    }
}


// ************************************************************************* //
