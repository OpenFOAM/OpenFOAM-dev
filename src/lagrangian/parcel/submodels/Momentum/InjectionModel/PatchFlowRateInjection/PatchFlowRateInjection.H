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

Class
    Foam::PatchFlowRateInjection

Description
    Patch injection, by using patch flow rate to determine concentration and
    velocity.

    User specifies:
      - Total mass to inject
      - Name of patch
      - Injection duration
      - Injection target concentration/carrier volume flow rate

    Properties:
      - Initial parcel velocity given by local flow velocity
      - Parcel diameters obtained by distribution model
      - Parcels injected randomly across the patch

SourceFiles
    PatchFlowRateInjection.C

\*---------------------------------------------------------------------------*/

#ifndef PatchFlowRateInjection_H
#define PatchFlowRateInjection_H

#include "InjectionModel.H"
#include "patchInjectionBase.H"
#include "TimeFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class distributionModel;

/*---------------------------------------------------------------------------*\
                   Class PatchFlowRateInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PatchFlowRateInjection
:
    public InjectionModel<CloudType>,
    public patchInjectionBase
{
    // Private Data

        //- Name of carrier (mass or volume) flux field
        const word phiName_;

        //- Name of carrier density field
        const word rhoName_;

        //- Injection duration [s]
        scalar duration_;

        //- Concentration profile of particle volume to carrier volume [-]
        const TimeFunction1<scalar> concentration_;

        //- Parcels to introduce per unit volume flow rate m3 [n/m^3]
        const scalar parcelConcentration_;

        //- Parcel size distribution model
        const autoPtr<distributionModel> sizeDistribution_;


public:

    //- Runtime type information
    TypeName("patchFlowRateInjection");


    // Constructors

        //- Construct from dictionary
        PatchFlowRateInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        PatchFlowRateInjection(const PatchFlowRateInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new PatchFlowRateInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~PatchFlowRateInjection();


    // Member Functions

        //- Inherit updateMesh from patchInjectionBase
        using patchInjectionBase::updateMesh;

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        scalar timeEnd() const;

        //- Return the total volumetric flow rate across the patch [m^3/s]
        virtual scalar flowRate() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Inherit setPositionAndCell from patchInjectionBase
            using patchInjectionBase::setPositionAndCell;

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFacei,
                label& tetPti
            );

            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PatchFlowRateInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
