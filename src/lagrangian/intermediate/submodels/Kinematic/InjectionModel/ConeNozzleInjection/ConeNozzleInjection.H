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

Class
    Foam::ConeNozzleInjection

Description
    Cone injection

    - User specifies
      - time of start of injection
      - injector position
      - direction (along injection axis)
      - parcel flow rate
      - inner and outer half-cone angles

    - Parcel diameters obtained by size distribution model

    - Parcel velocity is calculated as:

        - Constant velocity
            U = <specified by user>
        - Pressure driven velocity
            U = sqrt(2*(Pinj - Pamb)/rho)
        - Flow rate and discharge
            U = V_dot/(A*Cd)

SourceFiles
    ConeNozzleInjection.C

\*---------------------------------------------------------------------------*/

#ifndef ConeNozzleInjection_H
#define ConeNozzleInjection_H

#include "InjectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class Type>
class TimeDataEntry;

class distributionModel;

/*---------------------------------------------------------------------------*\
                     Class ConeNozzleInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ConeNozzleInjection
:
    public InjectionModel<CloudType>
{
public:

    //- Injection method enumeration
    enum injectionMethod
    {
        imPoint,
        imDisc
    };

    //- Flow type enumeration
    enum flowType
    {
        ftConstantVelocity,
        ftPressureDrivenVelocity,
        ftFlowRateAndDischarge
    };


private:

    // Private data

        //- Point/disc injection method
        injectionMethod injectionMethod_;

        //- Flow type
        flowType flowType_;

        //- Outer nozzle diameter [m]
        const scalar outerDiameter_;

        //- Inner nozzle diameter [m]
        const scalar innerDiameter_;

        //- Injection duration [s]
        scalar duration_;

        //- Injector position [m]
        vector position_;

        //- Cell containing injector position []
        label injectorCell_;

        //- Index of tet face for injector cell
        label tetFaceI_;

        //- Index of tet point for injector cell
        label tetPtI_;

        //- Injector direction []
        vector direction_;

        //- Number of parcels to introduce per second []
        const label parcelsPerSecond_;

        //- Flow rate profile relative to SOI []
        const TimeDataEntry<scalar> flowRateProfile_;

        //- Inner half-cone angle relative to SOI [deg]
        const TimeDataEntry<scalar> thetaInner_;

        //- Outer half-cone angle relative to SOI [deg]
        const TimeDataEntry<scalar> thetaOuter_;

        //- Parcel size PDF model
        const autoPtr<distributionModels::distributionModel> sizeDistribution_;


        // Tangential vectors to the direction vector

            //- First tangential vector
            vector tanVec1_;

            //- Second tangential vector
            vector tanVec2_;

            //- injection vector orthogonal to direction
            vector normal_;


        // Velocity model coefficients

            //- Constant velocity [m/s]
            scalar UMag_;

            //- Discharge coefficient, relative to SOI [m/s]
            TimeDataEntry<scalar> Cd_;

            //- Injection pressure [Pa]
            TimeDataEntry<scalar> Pinj_;


    // Private Member Functions

        //- Set the injection type
        void setInjectionMethod();

        //- Set the injection flow type
        void setFlowType();


public:

    //- Runtime type information
    TypeName("coneNozzleInjection");


    // Constructors

        //- Construct from dictionary
        ConeNozzleInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ConeNozzleInjection(const ConeNozzleInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType> > clone() const
        {
            return autoPtr<InjectionModel<CloudType> >
            (
                new ConeNozzleInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ConeNozzleInjection();


    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Set the injection position and owner cell
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFaceI,
                label& tetPtI
            );

            //- Set the parcel properties
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
#   include "ConeNozzleInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
