/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "SurfaceFilmModel.H"
#include "surfaceFilmModel.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::SurfaceFilmModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
    g_(owner.g()),
    ejectedParcelType_(0),
    massParcelPatch_(0),
    diameterParcelPatch_(0),
    UFilmPatch_(0),
    rhoFilmPatch_(0),
    deltaFilmPatch_(0),
    nParcelsTransferred_(0),
    nParcelsInjected_(0)
{}


template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::SurfaceFilmModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
    g_(owner.g()),
    ejectedParcelType_
    (
        this->coeffDict().lookupOrDefault("ejectedParcelType", -1)
    ),
    massParcelPatch_(0),
    diameterParcelPatch_(0),
    UFilmPatch_(0),
    rhoFilmPatch_(0),
    deltaFilmPatch_(owner.mesh().boundary().size()),
    nParcelsTransferred_(0),
    nParcelsInjected_(0)
{}


template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::SurfaceFilmModel
(
    const SurfaceFilmModel<CloudType>& sfm
)
:
    CloudSubModelBase<CloudType>(sfm),
    g_(sfm.g_),
    ejectedParcelType_(sfm.ejectedParcelType_),
    massParcelPatch_(sfm.massParcelPatch_),
    diameterParcelPatch_(sfm.diameterParcelPatch_),
    UFilmPatch_(sfm.UFilmPatch_),
    rhoFilmPatch_(sfm.rhoFilmPatch_),
    deltaFilmPatch_(sfm.deltaFilmPatch_),
    nParcelsTransferred_(sfm.nParcelsTransferred_),
    nParcelsInjected_(sfm.nParcelsInjected_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::~SurfaceFilmModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
template<class TrackData>
void Foam::SurfaceFilmModel<CloudType>::inject(TrackData& td)
{
    if (!this->active())
    {
        return;
    }

    // Retrieve the film model from the owner database
    const regionModels::surfaceFilmModels::surfaceFilmModel& filmModel =
        this->owner().mesh().time().objectRegistry::template lookupObject
        <regionModels::surfaceFilmModels::surfaceFilmModel>
        (
            "surfaceFilmProperties"
        );

    if (!filmModel.active())
    {
        return;
    }

    const labelList& filmPatches = filmModel.intCoupledPatchIDs();
    const labelList& primaryPatches = filmModel.primaryPatchIDs();

    const fvMesh& mesh = this->owner().mesh();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    forAll(filmPatches, i)
    {
        const label filmPatchi = filmPatches[i];
        const label primaryPatchi = primaryPatches[i];

        const labelList& injectorCellsPatch = pbm[primaryPatchi].faceCells();

        cacheFilmFields(filmPatchi, primaryPatchi, filmModel);

        const vectorField& Cf = mesh.C().boundaryField()[primaryPatchi];
        const vectorField& Sf = mesh.Sf().boundaryField()[primaryPatchi];
        const scalarField& magSf = mesh.magSf().boundaryField()[primaryPatchi];

        forAll(injectorCellsPatch, j)
        {
            if (diameterParcelPatch_[j] > 0)
            {
                const label celli = injectorCellsPatch[j];

                // The position could bein any tet of the decomposed cell,
                // so arbitrarily choose the first face of the cell as the
                // tetFace and the first point on the face after the base
                // point as the tetPt.  The tracking will pick the cell
                // consistent with the motion in the first tracking step.
                const label tetFacei = this->owner().mesh().cells()[celli][0];
                const label tetPtI = 1;

//                const point& pos = this->owner().mesh().C()[celli];

                const scalar offset =
                    max
                    (
                        diameterParcelPatch_[j],
                        deltaFilmPatch_[primaryPatchi][j]
                    );
                const point pos = Cf[j] - 1.1*offset*Sf[j]/magSf[j];

                // Create a new parcel
                parcelType* pPtr =
                    new parcelType
                    (
                        this->owner().pMesh(),
                        pos,
                        celli,
                        tetFacei,
                        tetPtI
                    );

                // Check/set new parcel thermo properties
                td.cloud().setParcelThermoProperties(*pPtr, 0.0);

                setParcelProperties(*pPtr, j);

                if (pPtr->nParticle() > 0.001)
                {
                    // Check new parcel properties
    //                td.cloud().checkParcelProperties(*pPtr, 0.0, true);
                    td.cloud().checkParcelProperties(*pPtr, 0.0, false);

                    // Add the new parcel to the cloud
                    td.cloud().addParticle(pPtr);

                    nParcelsInjected_++;
                }
                else
                {
                    // TODO: cache mass and re-distribute?
                    delete pPtr;
                }
            }
        }
    }
}


template<class CloudType>
void Foam::SurfaceFilmModel<CloudType>::cacheFilmFields
(
    const label filmPatchi,
    const label primaryPatchi,
    const regionModels::surfaceFilmModels::surfaceFilmModel& filmModel
)
{
    massParcelPatch_ = filmModel.cloudMassTrans().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, massParcelPatch_);

    diameterParcelPatch_ =
        filmModel.cloudDiameterTrans().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, diameterParcelPatch_, maxEqOp<scalar>());

    UFilmPatch_ = filmModel.Us().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, UFilmPatch_);

    rhoFilmPatch_ = filmModel.rho().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, rhoFilmPatch_);

    deltaFilmPatch_[primaryPatchi] =
        filmModel.delta().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, deltaFilmPatch_[primaryPatchi]);
}


template<class CloudType>
void Foam::SurfaceFilmModel<CloudType>::setParcelProperties
(
    parcelType& p,
    const label filmFacei
) const
{
    // Set parcel properties
    scalar vol = mathematical::pi/6.0*pow3(diameterParcelPatch_[filmFacei]);
    p.d() = diameterParcelPatch_[filmFacei];
    p.U() = UFilmPatch_[filmFacei];
    p.rho() = rhoFilmPatch_[filmFacei];

    p.nParticle() = massParcelPatch_[filmFacei]/p.rho()/vol;

    if (ejectedParcelType_ >= 0)
    {
        p.typeId() = ejectedParcelType_;
    }
}


template<class CloudType>
void Foam::SurfaceFilmModel<CloudType>::info(Ostream& os)
{
    label nTrans0 =
        this->template getModelProperty<label>("nParcelsTransferred");

    label nInject0 =
        this->template getModelProperty<label>("nParcelsInjected");

    label nTransTotal =
        nTrans0 + returnReduce(nParcelsTransferred_, sumOp<label>());

    label nInjectTotal =
        nInject0 + returnReduce(nParcelsInjected_, sumOp<label>());

    os  << "    Parcels absorbed into film      = " << nTransTotal << nl
        << "    New film detached parcels       = " << nInjectTotal << endl;

    if (this->writeTime())
    {
        this->setModelProperty("nParcelsTransferred", nTransTotal);
        this->setModelProperty("nParcelsInjected", nInjectTotal);
        nParcelsTransferred_ = 0;
        nParcelsInjected_ = 0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SurfaceFilmModelNew.C"

// ************************************************************************* //
