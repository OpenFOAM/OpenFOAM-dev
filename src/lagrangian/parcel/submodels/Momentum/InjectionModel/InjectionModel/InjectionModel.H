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

Class
    Foam::InjectionModel

Description
    Templated injection model class.

    The injection model nominally describes the parcel:
    - position
    - diameter
    - velocity
    In this case, the fullyDescribed() flag should be set to 0 (false). When
    the parcel is then added to the cloud, the remaining properties are
    populated using values supplied in the constant properties.

    If, however, all of a parcel's properties are described in the model, the
    fullDescribed() flag should be set to 1 (true).


SourceFiles
    InjectionModel.C
    InjectionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef InjectionModel_H
#define InjectionModel_H

#include "injectionModel.H"
#include "CloudSubModelBase.H"
#include "particle.H"
#include "TimeFunction1.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class InjectionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class InjectionModel
:
    public injectionModel,
    public CloudSubModelBase<CloudType>
{
public:

    //- Convenience typedef for parcelType
    typedef typename CloudType::parcelType parcelType;


protected:

    // Protected data

        // Global injection properties

            //- Start of injection [s]
            scalar SOI_;

            //- Total mass injected to date [kg]
            scalar massInjected_;


        // Counters

            //- Number of injections counter
            label nInjections_;

            //- Running counter of total number of parcels added
            label parcelsAddedTotal_;


        // Injection properties per Lagrangian time step

            //- Fixed nParticle to assign to parcels. Only valid if
            //  uniformParcelSize is nParticle.
            scalar nParticleFixed_;

            //- Size uniform to all parcels
            uniformParcelSize uniformParcelSize_;

            //- Continuous phase time at start of injection time step [s]
            scalar time0_;

            //- Time at start of injection time step [s]
            scalar timeStep0_;


    // Protected Member Functions

        //- Read the total mass value for instantaneous injections
        scalar readMassTotal
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Read the duration for continuous injections
        scalar readDuration
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Read the mass flow rate function for continuous injections
        TimeFunction1<scalar> readMassFlowRate
        (
            const dictionary& dict,
            CloudType& owner,
            const scalar duration
        );

        //- Read the number of parcels injected per second for continuous
        //  injections
        TimeFunction1<scalar> readParcelsPerSecond
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Get the index of this injector
        label index() const;

        //- Find the cell that contains the supplied position
        //  Will modify position slightly towards the owner cell centroid to
        //  ensure that it lies in a cell and not edge/face
        bool findCellAtPosition
        (
            const point& position,
            barycentric& coordinates,
            label& celli,
            label& tetFacei,
            label& tetPti,
            bool errorOnNotFound = true
        );

        //- Constrain a parcel's position appropriately to the geometric
        //  dimensions of the mesh
        void constrainPosition
        (
            typename CloudType::parcelType::trackingData& td,
            typename CloudType::parcelType& parcel
        );

        //- Return the sampling moment to be used by the size distribution
        label sizeSampleQ() const;

        //- Set number of particles to inject given parcel properties
        void setNumberOfParticles
        (
            PtrList<parcelType>& parcelPtrs,
            const scalar mass
        ) const;

        //- Post injection checks
        void postInjectCheck
        (
            const label parcelsAdded,
            const scalar massAdded
        );


public:

    //- Runtime type information
    TypeName("injectionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        InjectionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelType
        ),
        (dict, owner, modelType)
    );


    // Constructors

        //- Construct null from owner
        InjectionModel(CloudType& owner);

        //- Construct from dictionary
        InjectionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName,
            const word& modelType
        );

        //- Construct copy
        InjectionModel(const InjectionModel<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~InjectionModel();


    // Selectors

        //- Selector with lookup from dictionary
        static autoPtr<InjectionModel<CloudType>> New
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Selector with name and type
        static autoPtr<InjectionModel<CloudType>> New
        (
            const dictionary& dict,
            const word& modelName,
            const word& modelType,
            CloudType& owner
        );


    // Member Functions

        // Mapping

            //- Update mesh
            virtual void topoChange();


        // Global information

            //- Return the start-of-injection time
            inline scalar timeStart() const;

            //- Return mass of particles injected (cumulative)
            inline scalar massInjected() const;

            //- Return the end-of-injection time
            virtual scalar timeEnd() const = 0;

            //- Number of parcels to introduce relative to SOI
            virtual label nParcelsToInject
            (
                const scalar time0,
                const scalar time1
            ) = 0;

            //- Parcel mass to introduce relative to SOI
            virtual scalar massToInject
            (
                const scalar time0,
                const scalar time1
            ) = 0;

            //- Return the average injected parcel mass
            scalar averageParcelMass();


            // Counters

                //- Return the number of injections
                inline label nInjections() const;

                //- Return the total number parcels added
                inline label parcelsAddedTotal() const;


        // Per-injection event functions

            //- Main injection loop
            template<class TrackCloudType>
            void inject
            (
                TrackCloudType& cloud,
                typename CloudType::parcelType::trackingData& td
            );

            //- Main injection loop - steady-state
            template<class TrackCloudType>
            void injectSteadyState
            (
                TrackCloudType& cloud,
                typename CloudType::parcelType::trackingData& td
            );


        // Injection geometry

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
            ) = 0;

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                parcelType& parcel
            ) = 0;

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const = 0;


        // I-O

            //- Write injection info to stream
            virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeInjectionModel(CloudType)                                          \
                                                                               \
    typedef Foam::CloudType::momentumCloudType momentumCloudType;            \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::InjectionModel<momentumCloudType>,                              \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            InjectionModel<momentumCloudType>,                                \
            dictionary                                                         \
        );                                                                     \
    }


#define makeInjectionModelType(SS, CloudType)                                  \
                                                                               \
    typedef Foam::CloudType::momentumCloudType momentumCloudType;            \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<momentumCloudType>, 0);      \
                                                                               \
    Foam::InjectionModel<momentumCloudType>::                                 \
        adddictionaryConstructorToTable<Foam::SS<momentumCloudType>>          \
            add##SS##CloudType##momentumCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "InjectionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "InjectionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
