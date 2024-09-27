/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::PatchInjection

Description
    Patch injection.

    User specifies:
      - Total mass to inject
      - Name of patch
      - Injection duration
      - Initial parcel velocity
      - Injection volume flow rate

    Properties:
      - Parcel diameters obtained by distribution model
      - Parcels injected randomly across the patch

SourceFiles
    PatchInjection.C

\*---------------------------------------------------------------------------*/

#ifndef PatchInjection_H
#define PatchInjection_H

#include "InjectionModel.H"
#include "patchInjectionBase.H"
#include "Function1.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class distribution;

/*---------------------------------------------------------------------------*\
                       Class PatchInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PatchInjection
:
    public InjectionModel<CloudType>,
    public patchInjectionBase
{
    // Private Data

        //- Injection duration [s]
        const scalar duration_;

        //- Mass flow rate relative to SOI [kg/s]
        autoPtr<Function1<scalar>> massFlowRate_;

        //- Number of parcels to introduce per second [1/s]
        autoPtr<Function1<scalar>> parcelsPerSecond_;

        //- Fixed initial parcel velocity [m/s]
        vector U0_;

        //- Name of the field from which to interpolate the initial parcel
        //  velocity
        word U0Name_;

        //- Initial parcel interpolation engine. Invalid if using a fixed
        //  parcel velocity. Contains a null reference if using the carrier
        //  velocity (as the track data interpolation can be used). Contains an
        //  actual interpolation engine if using a different velocity.
        tmpNrc<interpolation<vector>> U0InterpPtr_;

        //- Parcel size distribution model
        const autoPtr<distribution> sizeDistribution_;


protected:

    // Protected Member Functions

        //- Pre injection hook
        virtual void preInject
        (
            typename CloudType::parcelType::trackingData& td
        );

        //- Post injection hook
        virtual void postInject
        (
            const label parcelsAdded,
            const scalar massAdded,
            typename CloudType::parcelType::trackingData& td
        );


public:

    //- Runtime type information
    TypeName("patchInjection");


    // Constructors

        //- Construct from dictionary
        PatchInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        PatchInjection(const PatchInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new PatchInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~PatchInjection();


    // Member Functions

        //- Inherit topoChange from patchInjectionBase
        using patchInjectionBase::topoChange;

        //- Set injector locations when mesh is updated
        virtual void topoChange();

        //- Return the end-of-injection time
        virtual scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual scalar nParcelsToInject(const scalar time0, const scalar time1);

        //- Parcel mass to introduce relative to SOI
        virtual scalar massToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Inherit setPositionAndCell from patchInjectionBase
            using patchInjectionBase::setPositionAndCell;

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                barycentric& coordinates,
                label& celli,
                label& tetFacei,
                label& tetPti,
                label& facei
            );

            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType::trackingData& td,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PatchInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
