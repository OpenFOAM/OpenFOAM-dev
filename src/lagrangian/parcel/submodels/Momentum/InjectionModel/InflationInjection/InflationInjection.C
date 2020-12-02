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

#include "InflationInjection.H"
#include "mathematicalConstants.H"
#include "PackedBoolList.H"
#include "cellSet.H"
#include "ListListOps.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InflationInjection<CloudType>::InflationInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    generationSetName_(this->coeffDict().lookup("generationCellSet")),
    inflationSetName_(this->coeffDict().lookup("inflationCellSet")),
    generationCells_(),
    inflationCells_(),
    duration_(this->coeffDict().template lookup<scalar>("duration")),
    flowRateProfile_
    (
        TimeFunction1<scalar>
        (
            owner.db().time(),
            "flowRateProfile",
            this->coeffDict()
        )
    ),
    growthRate_
    (
        TimeFunction1<scalar>
        (
            owner.db().time(),
            "growthRate",
            this->coeffDict()
        )
    ),
    newParticles_(),
    volumeAccumulator_(0.0),
    fraction_(1.0),
    selfSeed_(this->coeffDict().lookupOrDefault("selfSeed", false)),
    dSeed_(small),
    sizeDistribution_
    (
        distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"),
            owner.rndGen()
        )
    )
{
    duration_ = owner.db().time().userTimeToTime(duration_);

    if (selfSeed_)
    {
        dSeed_ = this->coeffDict().template lookup<scalar>("dSeed");
    }

    cellSet generationCells(this->owner().mesh(), generationSetName_);

    generationCells_ = generationCells.toc();

    cellSet inflationCells(this->owner().mesh(), inflationSetName_);

    // Union of cellSets
    inflationCells |= generationCells;

    inflationCells_ = inflationCells.toc();

    if (Pstream::parRun())
    {
        scalar generationVolume = 0.0;

        forAll(generationCells_, gCI)
        {
            label cI = generationCells_[gCI];

            generationVolume += this->owner().mesh().cellVolumes()[cI];
        }

        scalar totalGenerationVolume = generationVolume;

        reduce(totalGenerationVolume, sumOp<scalar>());

        fraction_ = generationVolume/totalGenerationVolume;
    }

    // Set total volume/mass to inject
    this->volumeTotal_ = fraction_*flowRateProfile_.integral(0.0, duration_);
    this->massTotal_ *= fraction_;
}


template<class CloudType>
Foam::InflationInjection<CloudType>::InflationInjection
(
    const Foam::InflationInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    generationSetName_(im.generationSetName_),
    inflationSetName_(im.inflationSetName_),
    generationCells_(im.generationCells_),
    inflationCells_(im.inflationCells_),
    duration_(im.duration_),
    flowRateProfile_(im.flowRateProfile_),
    growthRate_(im.growthRate_),
    newParticles_(im.newParticles_),
    volumeAccumulator_(im.volumeAccumulator_),
    fraction_(im.fraction_),
    selfSeed_(im.selfSeed_),
    dSeed_(im.dSeed_),
    sizeDistribution_(im.sizeDistribution_().clone().ptr())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InflationInjection<CloudType>::~InflationInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::InflationInjection<CloudType>::updateMesh()
{}


template<class CloudType>
Foam::scalar Foam::InflationInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::InflationInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    const polyMesh& mesh = this->owner().mesh();

    List<DynamicList<typename CloudType::parcelType*>>& cellOccupancy =
        this->owner().cellOccupancy();

    scalar gR = growthRate_.value(time1);

    scalar dT = time1 - time0;

    // Inflate existing particles

    forAll(inflationCells_, iCI)
    {
        label cI = inflationCells_[iCI];

        typename CloudType::parcelType* pPtr = nullptr;

        forAll(cellOccupancy[cI], cPI)
        {
            pPtr = cellOccupancy[cI][cPI];

            scalar dTarget = pPtr->dTarget();

            pPtr->d() = min(dTarget, pPtr->d() + gR*dT);
        }
    }

    // Generate new particles

    newParticles_.clear();

    Random& rnd = this->owner().rndGen();

    // Diameter factor, when splitting particles into 4, this is the
    // factor that modifies the diameter.
    scalar dFact = sqrt(2.0)/(sqrt(3.0) + sqrt(2.0));

    if ((time0 >= 0.0) && (time0 < duration_))
    {
        volumeAccumulator_ +=
            fraction_*flowRateProfile_.integral(time0, time1);
    }

    labelHashSet cellCentresUsed;

    // Loop escape counter
    label maxIterations = max
    (
        1,
        (10*volumeAccumulator_)
       /CloudType::parcelType::volume(sizeDistribution_().minValue())
    );

    label iterationNo = 0;

    // Info<< "Accumulated volume to inject: "
    //     << returnReduce(volumeAccumulator_, sumOp<scalar>()) << endl;

    while (!generationCells_.empty() && volumeAccumulator_ > 0)
    {
        if (iterationNo > maxIterations)
        {
            WarningInFunction
                << "Maximum particle split iterations ("
                << maxIterations << ") exceeded" << endl;

            break;
        }

        label cI =
            generationCells_[rnd.sampleAB<label>(0, generationCells_.size())];

        // Pick a particle at random from the cell - if there are
        // none, insert one at the cell centre.  Otherwise, split an
        // existing particle into four new ones.

        if (cellOccupancy[cI].empty())
        {
            if (selfSeed_ && !cellCentresUsed.found(cI))
            {
                scalar dNew = sizeDistribution_().sample();

                newParticles_.append
                (
                    vectorPairScalarPair
                    (
                        Pair<vector>(mesh.cellCentres()[cI], Zero),
                        Pair<scalar>(dSeed_, dNew)
                    )
                );

                volumeAccumulator_ -= CloudType::parcelType::volume(dNew);

                cellCentresUsed.insert(cI);
            }
        }
        else
        {
            label cPI = rnd.sampleAB<label>(0, cellOccupancy[cI].size());

            // This has to be a reference to the pointer so that it
            // can be set to nullptr when the particle is deleted.
            typename CloudType::parcelType*& pPtr = cellOccupancy[cI][cPI];

            if (pPtr != nullptr)
            {
                scalar pD = pPtr->d();

                // Select bigger particles by preference
                if ((pD/pPtr->dTarget()) < rnd.sample01<scalar>())
                {
                    continue;
                }

                const point& pP = pPtr->position();
                const vector& pU = pPtr->U();

                // Generate a tetrahedron of new positions with the
                // four new spheres fitting inside the old one, where
                // a is the diameter of the new spheres, and is
                // related to the diameter of the enclosing sphere, A,
                // by a = sqrt(2)*A/(sqrt(3) + sqrt(2));

                // Positions around the origin, which is the
                // tetrahedron centroid (centre of old sphere).

                // x = a/sqrt(3)
                // r = a/(2*sqrt(6))
                // R = sqrt(3)*a/(2*sqrt(2))
                // d = a/(2*sqrt(3))

                // p0(x, 0, -r)
                // p1(-d, a/2, -r)
                // p2(-d, -a/2, -r)
                // p3(0, 0, R)

                scalar a = pD*dFact;

                scalar x = a/sqrt(3.0);
                scalar r = a/(2.0*sqrt(6.0));
                scalar R = sqrt(3.0)*a/(2.0*sqrt(2.0));
                scalar d = a/(2.0*sqrt(3.0));

                scalar dNew = sizeDistribution_().sample();
                scalar volNew = CloudType::parcelType::volume(dNew);

                newParticles_.append
                (
                    vectorPairScalarPair
                    (
                        Pair<vector>(vector(x, 0, -r)  + pP, pU),
                        Pair<scalar>(a, dNew)
                    )
                );
                volumeAccumulator_ -= volNew;

                dNew = sizeDistribution_().sample();
                newParticles_.append
                (
                    vectorPairScalarPair
                    (
                        Pair<vector>(vector(-d, a/2, -r) + pP, pU),
                        Pair<scalar>(a, dNew)
                    )
                );
                volumeAccumulator_ -= volNew;

                dNew = sizeDistribution_().sample();
                newParticles_.append
                (
                    vectorPairScalarPair
                    (
                        Pair<vector>(vector(-d, -a/2, -r) + pP, pU),
                        Pair<scalar>(a, dNew)
                    )
                );
                volumeAccumulator_ -= volNew;

                dNew = sizeDistribution_().sample();
                newParticles_.append
                (
                    vectorPairScalarPair
                    (
                        Pair<vector>(vector(0, 0, R) + pP, pU),
                        Pair<scalar>(a, dNew)
                    )
                );
                volumeAccumulator_ -= volNew;

                // Account for the lost volume of the particle which
                // is to be deleted
                volumeAccumulator_ += CloudType::parcelType::volume
                (
                    pPtr->dTarget()
                );

                this->owner().deleteParticle(*pPtr);

                pPtr = nullptr;
            }
        }

        iterationNo++;
    }

    if (Pstream::parRun())
    {
        List<List<vectorPairScalarPair>> gatheredNewParticles
        (
            Pstream::nProcs()
        );

        gatheredNewParticles[Pstream::myProcNo()] = newParticles_;

        // Gather data onto master
        Pstream::gatherList(gatheredNewParticles);

        // Combine
        List<vectorPairScalarPair> combinedNewParticles
        (
            ListListOps::combine<List<vectorPairScalarPair>>
            (
                gatheredNewParticles,
                accessOp<List<vectorPairScalarPair>>()
            )
        );

        if (Pstream::master())
        {
            newParticles_ = combinedNewParticles;
        }

        Pstream::scatter(newParticles_);
    }

    return newParticles_.size();
}


template<class CloudType>
Foam::scalar Foam::InflationInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return fraction_*flowRateProfile_.integral(time0, time1);
    }
    else
    {
        return 0.0;
    }
}


template<class CloudType>
void Foam::InflationInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label,
    const scalar,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    position = newParticles_[parcelI].first().first();

    this->findCellAtPosition
    (
        cellOwner,
        tetFacei,
        tetPti,
        position,
        false
    );
}


template<class CloudType>
void Foam::InflationInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    parcel.U() = newParticles_[parcelI].first().second();

    parcel.d() = newParticles_[parcelI].second().first();

    parcel.dTarget() = newParticles_[parcelI].second().second();
}


template<class CloudType>
bool Foam::InflationInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::InflationInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
