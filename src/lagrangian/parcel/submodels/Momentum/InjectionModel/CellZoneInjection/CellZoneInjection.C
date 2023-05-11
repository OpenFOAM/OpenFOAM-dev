/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "CellZoneInjection.H"
#include "mathematicalConstants.H"
#include "polyMeshTetDecomposition.H"
#include "globalIndex.H"
#include "Pstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::CellZoneInjection<CloudType>::setPositions
(
    const labelList& cellZoneCells
)
{
    const fvMesh& mesh = this->owner().mesh();

    Random& rnd = this->owner().rndGen();

    const label nCells = cellZoneCells.size();
    DynamicList<barycentric> injectorCoordinates(nCells);
    DynamicList<label> injectorCells(nCells);
    DynamicList<label> injectorTetFaces(nCells);
    DynamicList<label> injectorTetPts(nCells);

    scalar newParticlesTotal = 0;
    label addParticlesTotal = 0;

    forAll(cellZoneCells, i)
    {
        const label celli = cellZoneCells[i];

        // Calc number of particles to add
        const scalar newParticles = mesh.V()[celli]*numberDensity_;
        newParticlesTotal += newParticles;
        label addParticles = floor(newParticles);
        addParticlesTotal += addParticles;

        const scalar diff = newParticlesTotal - addParticlesTotal;
        if (diff > 1)
        {
            label corr = floor(diff);
            addParticles += corr;
            addParticlesTotal += corr;
        }

        // Construct cell tet indices
        const List<tetIndices> cellTetIs =
            polyMeshTetDecomposition::cellTetIndices(mesh, celli);

        // Construct cell tet volume fractions
        scalarField cTetVFrac(cellTetIs.size());
        cTetVFrac[0] = cellTetIs[0].tet(mesh).mag();
        for (label tetI = 1; tetI < cellTetIs.size(); tetI++)
        {
            cTetVFrac[tetI] =
                cTetVFrac[tetI-1] + cellTetIs[tetI].tet(mesh).mag();
        }
        cTetVFrac /= cTetVFrac.last();

        // Set new particle position and cellId
        for (label pI = 0; pI < addParticles; pI++)
        {
            const scalar volFrac = rnd.sample01<scalar>();
            label tetI = 0;
            forAll(cTetVFrac, vfI)
            {
                if (cTetVFrac[vfI] > volFrac)
                {
                    tetI = vfI;
                    break;
                }
            }

            injectorCoordinates.append(barycentric01(rnd));
            injectorCells.append(celli);
            injectorTetFaces.append(cellTetIs[tetI].face());
            injectorTetPts.append(cellTetIs[tetI].tetPt());
        }
    }

    // Accumulate into global lists. Set the coordinates and topology for local
    // particles, and leave remote ones with invalid data.
    globalIndex globalPositions(injectorCoordinates.size());

    List<barycentric> allInjectorCoordinates
    (
        globalPositions.size(),
        barycentric::uniform(NaN)
    );
    List<label> allInjectorCells(globalPositions.size(), -1);
    List<label> allInjectorTetFaces(globalPositions.size(), -1);
    List<label> allInjectorTetPts(globalPositions.size(), -1);

    SubList<barycentric>
    (
        allInjectorCoordinates,
        globalPositions.localSize(Pstream::myProcNo()),
        globalPositions.offset(Pstream::myProcNo())
    ) = injectorCoordinates;
    SubList<label>
    (
        allInjectorCells,
        globalPositions.localSize(Pstream::myProcNo()),
        globalPositions.offset(Pstream::myProcNo())
    ) = injectorCells;
    SubList<label>
    (
        allInjectorTetFaces,
        globalPositions.localSize(Pstream::myProcNo()),
        globalPositions.offset(Pstream::myProcNo())
    ) = injectorTetFaces;
    SubList<label>
    (
        allInjectorTetPts,
        globalPositions.localSize(Pstream::myProcNo()),
        globalPositions.offset(Pstream::myProcNo())
    ) = injectorTetPts;

    // Transfer data
    injectorCoordinates_.transfer(allInjectorCoordinates);
    injectorCells_.transfer(allInjectorCells);
    injectorTetFaces_.transfer(allInjectorTetFaces);
    injectorTetPts_.transfer(allInjectorTetPts);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CellZoneInjection<CloudType>::CellZoneInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    cellZoneName_(this->coeffDict().lookup("cellZone")),
    massTotal_(this->readMassTotal(dict, owner)),
    numberDensity_(this->coeffDict().template lookup<scalar>("numberDensity")),
    injectorCoordinates_(),
    injectorCells_(),
    injectorTetFaces_(),
    injectorTetPts_(),
    diameters_(),
    U0_(this->coeffDict().lookup("U0")),
    sizeDistribution_
    (
        distribution::New
        (
            this->coeffDict().subDict("sizeDistribution"),
            owner.rndGen(),
            this->sizeSampleQ()
        )
    )
{
    topoChange();
}


template<class CloudType>
Foam::CellZoneInjection<CloudType>::CellZoneInjection
(
    const CellZoneInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    cellZoneName_(im.cellZoneName_),
    massTotal_(im.massTotal_),
    numberDensity_(im.numberDensity_),
    injectorCoordinates_(im.injectorCoordinates_),
    injectorCells_(im.injectorCells_),
    injectorTetFaces_(im.injectorTetFaces_),
    injectorTetPts_(im.injectorTetPts_),
    diameters_(im.diameters_),
    U0_(im.U0_),
    sizeDistribution_(im.sizeDistribution_().clone().ptr())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CellZoneInjection<CloudType>::~CellZoneInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CellZoneInjection<CloudType>::topoChange()
{
    // Set/cache the injector cells
    const fvMesh& mesh = this->owner().mesh();
    const label zoneI = mesh.cellZones().findZoneID(cellZoneName_);

    if (zoneI < 0)
    {
        FatalErrorInFunction
            << "Unknown cell zone name: " << cellZoneName_
            << ". Valid cell zones are: " << mesh.cellZones().names()
            << nl << exit(FatalError);
    }

    const labelList& cellZoneCells = mesh.cellZones()[zoneI];
    const label nCells = cellZoneCells.size();
    const scalar nCellsTotal = returnReduce(nCells, sumOp<label>());
    const scalar VCells = sum(scalarField(mesh.V(), cellZoneCells));
    const scalar VCellsTotal = returnReduce(VCells, sumOp<scalar>());
    Info<< "    cell zone size      = " << nCellsTotal << endl;
    Info<< "    cell zone volume    = " << VCellsTotal << endl;

    if ((nCellsTotal == 0) || (VCellsTotal*numberDensity_ < 1))
    {
        WarningInFunction
            << "Number of particles to be added to cellZone " << cellZoneName_
            << " is zero" << endl;
    }
    else
    {
        setPositions(cellZoneCells);

        Info<< "    number density      = " << numberDensity_ << nl
            << "    number of particles = " << injectorCoordinates_.size()
            << endl;

        // Construct parcel diameters
        diameters_.setSize(injectorCoordinates_.size());
        forAll(diameters_, i)
        {
            diameters_[i] = sizeDistribution_->sample();
        }
    }
}


template<class CloudType>
Foam::scalar Foam::CellZoneInjection<CloudType>::timeEnd() const
{
    // Not used
    return this->SOI_;
}


template<class CloudType>
Foam::label Foam::CellZoneInjection<CloudType>::nParcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    // All parcels introduced at SOI
    if (0 >= time0 && 0 < time1)
    {
        return injectorCoordinates_.size();
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::CellZoneInjection<CloudType>::massToInject
(
    const scalar time0,
    const scalar time1
)
{
    // All parcels introduced at SOI
    if (0 >= time0 && 0 < time1)
    {
        return massTotal_;
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
void Foam::CellZoneInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label nParcels,
    const scalar time,
    barycentric& coordinates,
    label& celli,
    label& tetFacei,
    label& tetPti,
    label& facei
)
{
    coordinates = injectorCoordinates_[parcelI];
    celli = injectorCells_[parcelI];
    tetFacei = injectorTetFaces_[parcelI];
    tetPti = injectorTetPts_[parcelI];
}


template<class CloudType>
void Foam::CellZoneInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // set particle velocity
    parcel.U() = U0_;

    // set particle diameter
    parcel.d() = diameters_[parcelI];
}


template<class CloudType>
bool Foam::CellZoneInjection<CloudType>::fullyDescribed() const
{
    return false;
}


// ************************************************************************* //
