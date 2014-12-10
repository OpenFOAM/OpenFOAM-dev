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

#include "InjectionModel.H"
#include "mathematicalConstants.H"
#include "meshTools.H"
#include "volFields.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class CloudType>
bool Foam::InjectionModel<CloudType>::validInjection(const label parcelI)
{
    notImplemented
    (
        "bool Foam::InjectionModel<CloudType>::validInjection(const label)"
    );

    return false;
}


template<class CloudType>
bool Foam::InjectionModel<CloudType>::prepareForNextTimeStep
(
    const scalar time,
    label& newParcels,
    scalar& newVolumeFraction
)
{
    // Initialise values
    newParcels = 0;
    newVolumeFraction = 0.0;
    bool validInjection = false;

    // Return if not started injection event
    if (time < SOI_)
    {
        timeStep0_ = time;
        return validInjection;
    }

    // Make times relative to SOI
    scalar t0 = timeStep0_ - SOI_;
    scalar t1 = time - SOI_;

    // Number of parcels to inject
    newParcels = this->parcelsToInject(t0, t1);

    // Volume of parcels to inject
    newVolumeFraction =
        this->volumeToInject(t0, t1)
       /(volumeTotal_ + ROOTVSMALL);

    if (newVolumeFraction > 0)
    {
        if (newParcels > 0)
        {
            timeStep0_ = time;
            validInjection = true;
        }
        else
        {
            // injection should have started, but not sufficient volume to
            // produce (at least) 1 parcel - hold value of timeStep0_
            validInjection = false;
        }
    }
    else
    {
        timeStep0_ = time;
        validInjection = false;
    }

    return validInjection;
}


template<class CloudType>
bool Foam::InjectionModel<CloudType>::findCellAtPosition
(
    label& cellI,
    label& tetFaceI,
    label& tetPtI,
    vector& position,
    bool errorOnNotFound
)
{
    const volVectorField& cellCentres = this->owner().mesh().C();

    const vector p0 = position;

    this->owner().mesh().findCellFacePt
    (
        position,
        cellI,
        tetFaceI,
        tetPtI
    );

    label procI = -1;

    if (cellI >= 0)
    {
        procI = Pstream::myProcNo();
    }

    reduce(procI, maxOp<label>());

    // Ensure that only one processor attempts to insert this Parcel

    if (procI != Pstream::myProcNo())
    {
        cellI = -1;
        tetFaceI = -1;
        tetPtI = -1;
    }

    // Last chance - find nearest cell and try that one - the point is
    // probably on an edge
    if (procI == -1)
    {
        cellI = this->owner().mesh().findNearestCell(position);

        if (cellI >= 0)
        {
            position += SMALL*(cellCentres[cellI] - position);

            if (this->owner().mesh().pointInCell(position, cellI))
            {
                procI = Pstream::myProcNo();
            }
        }

        reduce(procI, maxOp<label>());

        if (procI != Pstream::myProcNo())
        {
            cellI = -1;
            tetFaceI = -1;
            tetPtI = -1;
        }
    }

    if (procI == -1)
    {
        if (errorOnNotFound)
        {
            FatalErrorIn
            (
                "Foam::InjectionModel<CloudType>::findCellAtPosition"
                "("
                    "label&, "
                    "label&, "
                    "label&, "
                    "vector&, "
                    "bool"
                ")"
            )   << "Cannot find parcel injection cell. "
                << "Parcel position = " << p0 << nl
                << abort(FatalError);
        }
        else
        {
            return false;
        }
    }

    return true;
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::setNumberOfParticles
(
    const label parcels,
    const scalar volumeFraction,
    const scalar diameter,
    const scalar rho
)
{
    scalar nP = 0.0;
    switch (parcelBasis_)
    {
        case pbMass:
        {
            scalar volumep = pi/6.0*pow3(diameter);
            scalar volumeTot = massTotal_/rho;

            nP = (volumeFraction*volumeTot + delayedVolume_)/(parcels*volumep);
            break;
        }
        case pbNumber:
        {
            nP = massTotal_/(rho*volumeTotal_);
            break;
        }
        case pbFixed:
        {
            nP = nParticleFixed_;
            break;
        }
        default:
        {
            nP = 0.0;
            FatalErrorIn
            (
                "Foam::scalar "
                "Foam::InjectionModel<CloudType>::setNumberOfParticles"
                "("
                    "const label, "
                    "const scalar, "
                    "const scalar, "
                    "const scalar"
                ")"
            )<< "Unknown parcelBasis type" << nl
             << exit(FatalError);
        }
    }

    return nP;
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::postInjectCheck
(
    const label parcelsAdded,
    const scalar massAdded
)
{
    const label allParcelsAdded = returnReduce(parcelsAdded, sumOp<label>());

    if (allParcelsAdded > 0)
    {
        Info<< nl
            << "--> Cloud: " << this->owner().name()
            << " injector: " << this->modelName() << nl
            << "    Added " << allParcelsAdded << " new parcels" << nl << endl;
    }

    // Increment total number of parcels added
    parcelsAddedTotal_ += allParcelsAdded;

    // Increment total mass injected
    massInjected_ += returnReduce(massAdded, sumOp<scalar>());

    // Update time for start of next injection
    time0_ = this->owner().db().time().value();

    // Increment number of injections
    nInjections_++;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InjectionModel<CloudType>::InjectionModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
    SOI_(0.0),
    volumeTotal_(0.0),
    massTotal_(0.0),
    massFlowRate_(owner.db().time(), "massFlowRate"),
    massInjected_(this->template getModelProperty<scalar>("massInjected")),
    nInjections_(this->template getModelProperty<label>("nInjections")),
    parcelsAddedTotal_
    (
        this->template getModelProperty<scalar>("parcelsAddedTotal")
    ),
    parcelBasis_(pbNumber),
    nParticleFixed_(0.0),
    time0_(0.0),
    timeStep0_(this->template getModelProperty<scalar>("timeStep0")),
    delayedVolume_(0.0)
{}


template<class CloudType>
Foam::InjectionModel<CloudType>::InjectionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName,
    const word& modelType
)
:
    CloudSubModelBase<CloudType>(modelName, owner, dict, typeName, modelType),
    SOI_(0.0),
    volumeTotal_(0.0),
    massTotal_(0.0),
    massFlowRate_(owner.db().time(), "massFlowRate"),
    massInjected_(this->template getModelProperty<scalar>("massInjected")),
    nInjections_(this->template getModelProperty<scalar>("nInjections")),
    parcelsAddedTotal_
    (
        this->template getModelProperty<scalar>("parcelsAddedTotal")
    ),
    parcelBasis_(pbNumber),
    nParticleFixed_(0.0),
    time0_(owner.db().time().value()),
    timeStep0_(this->template getModelProperty<scalar>("timeStep0")),
    delayedVolume_(0.0)
{
    // Provide some info
    // - also serves to initialise mesh dimensions - needed for parallel runs
    //   due to lazy evaluation of valid mesh dimensions
    Info<< "    Constructing " << owner.mesh().nGeometricD() << "-D injection"
        << endl;

    if (owner.solution().transient())
    {
        this->coeffDict().lookup("massTotal") >> massTotal_;
        this->coeffDict().lookup("SOI") >> SOI_;
        SOI_ = owner.db().time().userTimeToTime(SOI_);
    }
    else
    {
        massFlowRate_.reset(this->coeffDict());
        massTotal_ = massFlowRate_.value(owner.db().time().value());
    }

    const word parcelBasisType = this->coeffDict().lookup("parcelBasisType");

    if (parcelBasisType == "mass")
    {
        parcelBasis_ = pbMass;
    }
    else if (parcelBasisType == "number")
    {
        parcelBasis_ = pbNumber;
    }
    else if (parcelBasisType == "fixed")
    {
        parcelBasis_ = pbFixed;

        Info<< "    Choosing nParticle to be a fixed value, massTotal "
            << "variable now does not determine anything."
            << endl;

        nParticleFixed_ = readScalar(this->coeffDict().lookup("nParticle"));
    }
    else
    {
        FatalErrorIn
        (
            "Foam::InjectionModel<CloudType>::InjectionModel"
            "("
                "const dictionary&, "
                "CloudType&, "
                "const word&"
            ")"
        )<< "parcelBasisType must be either 'number', 'mass' or 'fixed'" << nl
         << exit(FatalError);
    }
}


template<class CloudType>
Foam::InjectionModel<CloudType>::InjectionModel
(
    const InjectionModel<CloudType>& im
)
:
    CloudSubModelBase<CloudType>(im),
    SOI_(im.SOI_),
    volumeTotal_(im.volumeTotal_),
    massTotal_(im.massTotal_),
    massFlowRate_(im.massFlowRate_),
    massInjected_(im.massInjected_),
    nInjections_(im.nInjections_),
    parcelsAddedTotal_(im.parcelsAddedTotal_),
    parcelBasis_(im.parcelBasis_),
    nParticleFixed_(im.nParticleFixed_),
    time0_(im.time0_),
    timeStep0_(im.timeStep0_),
    delayedVolume_(im.delayedVolume_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InjectionModel<CloudType>::~InjectionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::InjectionModel<CloudType>::updateMesh()
{
    // do nothing
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::timeEnd() const
{
    notImplemented
    (
        "Foam::scalar Foam::InjectionModel<CloudType>::timeEnd() const"
    );

    return 0.0;
}


template<class CloudType>
Foam::label Foam::InjectionModel<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    notImplemented
    (
        "Foam::label Foam::InjectionModel<CloudType>::parcelsToInject"
        "("
            "const scalar, "
            "const scalar"
        ")"
    );

    return 0;
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    notImplemented
    (
        "Foam::scalar Foam::InjectionModel<CloudType>::volumeToInject"
        "("
            "const scalar, "
            "const scalar"
        ")"
    );

    return 0.0;
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::averageParcelMass()
{
    label nTotal = 0.0;
    if (this->owner().solution().transient())
    {
        nTotal = parcelsToInject(0.0, timeEnd() - timeStart());
    }
    else
    {
        nTotal = parcelsToInject(0.0, 1.0);
    }

    return massTotal_/nTotal;
}


template<class CloudType>
template<class TrackData>
void Foam::InjectionModel<CloudType>::inject(TrackData& td)
{
    if (!this->active())
    {
        return;
    }

    const scalar time = this->owner().db().time().value();

    // Prepare for next time step
    label parcelsAdded = 0;
    scalar massAdded = 0.0;
    label newParcels = 0;
    scalar newVolumeFraction = 0.0;

    if (prepareForNextTimeStep(time, newParcels, newVolumeFraction))
    {
        scalar delayedVolume = 0;

        const scalar trackTime = this->owner().solution().trackTime();
        const polyMesh& mesh = this->owner().mesh();
        typename TrackData::cloudType& cloud = td.cloud();

        // Duration of injection period during this timestep
        const scalar deltaT =
            max(0.0, min(trackTime, min(time - SOI_, timeEnd() - time0_)));

        // Pad injection time if injection starts during this timestep
        const scalar padTime = max(0.0, SOI_ - time0_);

        // Introduce new parcels linearly across carrier phase timestep
        for (label parcelI = 0; parcelI < newParcels; parcelI++)
        {
            if (validInjection(parcelI))
            {
                // Calculate the pseudo time of injection for parcel 'parcelI'
                scalar timeInj = time0_ + padTime + deltaT*parcelI/newParcels;

                // Determine the injection position and owner cell,
                // tetFace and tetPt
                label cellI = -1;
                label tetFaceI = -1;
                label tetPtI = -1;

                vector pos = vector::zero;

                setPositionAndCell
                (
                    parcelI,
                    newParcels,
                    timeInj,
                    pos,
                    cellI,
                    tetFaceI,
                    tetPtI
                );

                if (cellI > -1)
                {
                    // Lagrangian timestep
                    const scalar dt = time - timeInj;

                    // Apply corrections to position for 2-D cases
                    meshTools::constrainToMeshCentre(mesh, pos);

                    // Create a new parcel
                    parcelType* pPtr =
                        new parcelType(mesh, pos, cellI, tetFaceI, tetPtI);

                    // Check/set new parcel thermo properties
                    cloud.setParcelThermoProperties(*pPtr, dt);

                    // Assign new parcel properties in injection model
                    setProperties(parcelI, newParcels, timeInj, *pPtr);

                    // Check/set new parcel injection properties
                    cloud.checkParcelProperties(*pPtr, dt, fullyDescribed());

                    // Apply correction to velocity for 2-D cases
                    meshTools::constrainDirection
                    (
                        mesh,
                        mesh.solutionD(),
                        pPtr->U()
                    );

                    // Number of particles per parcel
                    pPtr->nParticle() =
                        setNumberOfParticles
                        (
                            newParcels,
                            newVolumeFraction,
                            pPtr->d(),
                            pPtr->rho()
                        );

                    if (pPtr->nParticle() >= 1.0)
                    {
                        parcelsAdded++;
                        massAdded += pPtr->nParticle()*pPtr->mass();

                        if (pPtr->move(td, dt))
                        {
                            td.cloud().addParticle(pPtr);
                        }
                        else
                        {
                            delete pPtr;
                        }
                    }
                    else
                    {
                        delayedVolume += pPtr->nParticle()*pPtr->volume();
                        delete pPtr;
                    }
                }
            }
        }

        delayedVolume_ = delayedVolume;
    }

    postInjectCheck(parcelsAdded, massAdded);
}


template<class CloudType>
template<class TrackData>
void Foam::InjectionModel<CloudType>::injectSteadyState
(
    TrackData& td,
    const scalar trackTime
)
{
    if (!this->active())
    {
        return;
    }

    const polyMesh& mesh = this->owner().mesh();
    typename TrackData::cloudType& cloud = td.cloud();

    massTotal_ = massFlowRate_.value(mesh.time().value());

    // Reset counters
    time0_ = 0.0;
    label parcelsAdded = 0;
    scalar massAdded = 0.0;

    // Set number of new parcels to inject based on first second of injection
    label newParcels = parcelsToInject(0.0, 1.0);

    // Inject new parcels
    for (label parcelI = 0; parcelI < newParcels; parcelI++)
    {
        // Volume to inject is split equally amongst all parcel streams
        scalar newVolumeFraction = 1.0/scalar(newParcels);

        // Determine the injection position and owner cell,
        // tetFace and tetPt
        label cellI = -1;
        label tetFaceI = -1;
        label tetPtI = -1;

        vector pos = vector::zero;

        setPositionAndCell
        (
            parcelI,
            newParcels,
            0.0,
            pos,
            cellI,
            tetFaceI,
            tetPtI
        );

        if (cellI > -1)
        {
            // Apply corrections to position for 2-D cases
            meshTools::constrainToMeshCentre(mesh, pos);

            // Create a new parcel
            parcelType* pPtr =
                new parcelType(mesh, pos, cellI, tetFaceI, tetPtI);

            // Check/set new parcel thermo properties
            cloud.setParcelThermoProperties(*pPtr, 0.0);

            // Assign new parcel properties in injection model
            setProperties(parcelI, newParcels, 0.0, *pPtr);

            // Check/set new parcel injection properties
            cloud.checkParcelProperties(*pPtr, 0.0, fullyDescribed());

            // Apply correction to velocity for 2-D cases
            meshTools::constrainDirection(mesh, mesh.solutionD(), pPtr->U());

            // Number of particles per parcel
            pPtr->nParticle() =
                setNumberOfParticles
                (
                    1,
                    newVolumeFraction,
                    pPtr->d(),
                    pPtr->rho()
                );

            // Add the new parcel
            td.cloud().addParticle(pPtr);

            massAdded += pPtr->nParticle()*pPtr->mass();
            parcelsAdded++;
        }
    }

    postInjectCheck(parcelsAdded, massAdded);
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label nParcels,
    const scalar time,
    vector& position,
    label& cellOwner,
    label& tetFaceI,
    label& tetPtI
)
{
    notImplemented
    (
        "void Foam::InjectionModel<CloudType>::setPositionAndCell"
        "("
            "const label, "
            "const label, "
            "const scalar, "
            "vector&, "
            "label&, "
            "label&, "
            "label&"
        ")"
    );
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::setProperties
(
    const label parcelI,
    const label nParcels,
    const scalar time,
    typename CloudType::parcelType& parcel
)
{
    notImplemented
    (
        "void Foam::InjectionModel<CloudType>::setProperties"
        "("
            "const label, "
            "const label, "
            "const scalar, "
            "typename CloudType::parcelType&"
        ")"
    );
}


template<class CloudType>
bool Foam::InjectionModel<CloudType>::fullyDescribed() const
{
    notImplemented
    (
        "bool Foam::InjectionModel<CloudType>::fullyDescribed() const"
    );

    return false;
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::info(Ostream& os)
{
    os  << "    " << this->modelName() << ":" << nl
        << "        number of parcels added     = " << parcelsAddedTotal_ << nl
        << "        mass introduced             = " << massInjected_ << nl;

    if (this->outputTime())
    {
        this->setModelProperty("massInjected", massInjected_);
        this->setModelProperty("nInjections", nInjections_);
        this->setModelProperty("parcelsAddedTotal", parcelsAddedTotal_);
        this->setModelProperty("timeStep0", timeStep0_);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "InjectionModelNew.C"

// ************************************************************************* //
