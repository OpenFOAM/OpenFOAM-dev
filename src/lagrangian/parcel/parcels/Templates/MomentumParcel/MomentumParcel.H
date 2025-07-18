/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::MomentumParcel

Description
    Momentum parcel class with rotational motion (as spherical
    particles only) and one/two-way coupling with the continuous
    phase.

    Sub-models include:
    - drag
    - turbulent dispersion
    - wall interactions

SourceFiles
    MomentumParcelI.H
    MomentumParcel.C
    MomentumParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef MomentumParcel_H
#define MomentumParcel_H

#include "particle.H"
#include "interpolation.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class MomentumParcel;

// Forward declaration of friend functions

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const MomentumParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                     Class MomentumParcelName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(MomentumParcel);


/*---------------------------------------------------------------------------*\
                       Class MomentumParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class MomentumParcel
:
    public ParcelType,
    public MomentumParcelName
{
    // Private Data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;

        //- Number of particle tracking attempts before we assume that it stalls
        static label maxTrackAttempts;


public:

    //- Class to hold momentum parcel constant properties
    class constantProperties
    {
    protected:

        // Protected data

            //- Constant properties dictionary
            const dictionary dict_;


    private:

        // Private Data

            //- Parcel type id - used for post-processing to flag the type
            //  of parcels issued by this cloud
            demandDrivenEntry<label> parcelTypeId_;

            //- Minimum density [kg/m^3]
            demandDrivenEntry<scalar> rhoMin_;

            //- Particle density [kg/m^3] (constant)
            demandDrivenEntry<scalar> rho0_;

            //- Minimum parcel mass [kg]
            demandDrivenEntry<scalar> minParcelMass_;


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);


        // Member Functions

            //- Return const access to the constant properties dictionary
            inline const dictionary& dict() const;

            //- Return const access to the parcel type id
            inline label parcelTypeId() const;

            //- Return const access to the minimum density
            inline scalar rhoMin() const;

            //- Return const access to the particle density
            inline scalar rho0() const;

            //- Return const access to the minimum parcel mass
            inline scalar minParcelMass() const;
    };


    class trackingData
    :
        public ParcelType::trackingData
    {
    private:

        // Private Data

            // Interpolators for continuous phase fields

                //- Density interpolator
                autoPtr<interpolation<scalar>> rhoInterp_;

                //- Velocity interpolator
                autoPtr<interpolation<vector>> UInterp_;

                //- Dynamic viscosity interpolator
                autoPtr<interpolation<scalar>> muInterp_;


            // Cached continuous phase properties

                //- Density [kg/m^3]
                scalar rhoc_;

                //- Velocity [m/s]
                vector Uc_;

                //- Viscosity [Pa.s]
                scalar muc_;


            //- Local gravitational or other body-force acceleration
            const vector& g_;


            // Tracking Parameters

                //- Track time. Total across the entire step.
                scalar trackTime_;

                //- Step fraction range to track between. Allows for the
                //  creation of sub-steps.
                Pair<scalar> stepFractionRange_;


    public:

        // Constructors

            //- Construct from components
            template <class TrackCloudType>
            inline trackingData(const TrackCloudType& cloud);


        // Member Functions

            //- Return const access to the interpolator for continuous
            //  phase density field
            inline const interpolation<scalar>& rhoInterp() const;

            //- Return const access to the interpolator for continuous
            //  phase velocity field
            inline const interpolation<vector>& UInterp() const;

            //- Return const access to the interpolator for continuous
            //  phase dynamic viscosity field
            inline const interpolation<scalar>& muInterp() const;

            //- Return the continuous phase density
            inline scalar rhoc() const;

            //- Access the continuous phase density
            inline scalar& rhoc();

            //- Return the continuous phase velocity
            inline const vector& Uc() const;

            //- Access the continuous phase velocity
            inline vector& Uc();

            //- Return the continuous phase viscosity
            inline scalar muc() const;

            //- Access the continuous phase viscosity
            inline scalar& muc();

            // Return the gravitational acceleration vector
            inline const vector& g() const;

            //- Return the tracking time
            inline scalar trackTime() const;

            //- Access the tracking time
            inline scalar& trackTime();

            //- Return the step fraction range to track between
            inline Pair<scalar> stepFractionRange() const;

            //- Access the step fraction range to track between
            inline Pair<scalar>& stepFractionRange();
    };


protected:

    // Protected data

        // Parcel properties

            //- Moving flag - tracking stopped when moving = false
            bool moving_;

            //- Parcel type id
            label typeId_;

            //- Number of particles in Parcel
            scalar nParticle_;

            //- Diameter [m]
            scalar d_;

            //- Target diameter [m]
            scalar dTarget_;

            //- Velocity of Parcel [m/s]
            vector U_;

            //- Density [kg/m^3]
            scalar rho_;

            //- Age [s]
            scalar age_;

            //- Time spent in turbulent eddy [s]
            scalar tTurb_;

            //- Turbulent velocity fluctuation [m/s]
            vector UTurb_;


    // Protected Member Functions

        //- Set cell values
        template<class TrackCloudType>
        void setCellValues(TrackCloudType& cloud, trackingData& td);

        //- Apply dispersion to the carrier phase velocity and update
        //  parcel turbulence parameters
        template<class TrackCloudType>
        void calcDispersion
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt
        );

        //- Correct cell values using latest transfer information
        template<class TrackCloudType>
        void cellValueSourceCorrection
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt
        );

        //- Update parcel properties over the time interval
        template<class TrackCloudType>
        void calc
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt
        );

        //- Calculate new particle velocity
        template<class TrackCloudType>
        const vector calcVelocity
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt,           // timestep
            const scalar Re,           // Reynolds number
            const scalar mu,           // local carrier viscosity
            const scalar mass,         // mass
            const vector& Su,          // explicit particle momentum source
            vector& dUTrans,           // momentum transfer to carrier
            scalar& Spu                // linearised drag coefficient
        ) const;


public:

    // Static Data Members

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " moving"
          + " typeId"
          + " nParticle"
          + " d"
          + " dTarget "
          + " (Ux Uy Uz)"
          + " rho"
          + " age"
          + " tTurb"
          + " (UTurbx UTurby UTurbz)"
        );


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline MomentumParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const label facei
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline MomentumParcel
        (
            const meshSearch& searchEngine,
            const vector& position,
            const label celli,
            label& nLocateBoundaryHits
        );

        //- Construct from Istream
        MomentumParcel(Istream& is, bool readFields = true);

        //- Construct as a copy
        MomentumParcel(const MomentumParcel& p);

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new MomentumParcel(*this));
        }

        //- Construct from Istream and return
        static autoPtr<MomentumParcel> New(Istream& is)
        {
            return autoPtr<MomentumParcel>(new MomentumParcel(is));
        }


    // Member Functions

        // Access

            //- Return const access to moving flag
            inline bool moving() const;

            //- Return const access to type id
            inline label typeId() const;

            //- Return const access to number of particles
            inline scalar nParticle() const;

            //- Return const access to diameter
            inline scalar d() const;

            //- Return const access to target diameter
            inline scalar dTarget() const;

            //- Return const access to velocity
            inline const vector& U() const;

            //- Return const access to density
            inline scalar rho() const;

            //- Return const access to the age
            inline scalar age() const;

            //- Return const access to time spent in turbulent eddy
            inline scalar tTurb() const;

            //- Return const access to turbulent velocity fluctuation
            inline const vector& UTurb() const;


        // Edit

            //- Return const access to moving flag
            inline bool& moving();

            //- Return access to type id
            inline label& typeId();

            //- Return access to number of particles
            inline scalar& nParticle();

            //- Return access to diameter
            inline scalar& d();

            //- Return access to target diameter
            inline scalar& dTarget();

            //- Return access to velocity
            inline vector& U();

            //- Return access to density
            inline scalar& rho();

            //- Return access to the age
            inline scalar& age();

            //- Return access to time spent in turbulent eddy
            inline scalar& tTurb();

            //- Return access to turbulent velocity fluctuation
            inline vector& UTurb();


        // Helper functions

            //- Cell owner mass
            inline scalar massCell(const trackingData& td) const;

            //- Particle mass
            inline scalar mass() const;

            //- Particle moment of inertia around diameter axis
            inline scalar momentOfInertia() const;

            //- Particle volume
            inline scalar volume() const;

            //- Particle volume for a given diameter
            inline static scalar volume(const scalar d);

            //- Particle projected area
            inline scalar areaP() const;

            //- Projected area for given diameter
            inline static scalar areaP(const scalar d);

            //- Particle surface area
            inline scalar areaS() const;

            //- Surface area for given diameter
            inline static scalar areaS(const scalar d);

            //- Reynolds number
            inline scalar Re(const trackingData& td) const;

            //- Reynolds number for given conditions
            inline static scalar Re
            (
                const scalar rhoc,
                const vector& U,
                const vector& Uc,
                const scalar d,
                const scalar muc
            );

            //- Weber number
            inline scalar We
            (
                const trackingData& td,
                const scalar sigma
            ) const;

            //- Weber number for given conditions
            inline static scalar We
            (
                const scalar rhoc,
                const vector& U,
                const vector& Uc,
                const scalar d,
                const scalar sigma
            );

            //- Eotvos number
            inline scalar Eo
            (
                const trackingData& td,
                const scalar sigma
            ) const;

            //- Eotvos number for given conditions
            inline static scalar Eo
            (
                const vector& g,
                const scalar rho,
                const scalar rhoc,
                const vector& U,
                const scalar d,
                const scalar sigma
            );


        // Tracking

            //- Move the parcel
            template<class TrackCloudType>
            bool move(TrackCloudType& cloud, trackingData& td);


        // Transformations

            //- Transform the physical properties of the particle
            //  according to the given transformation
            virtual void transformProperties(const transformer&);


        // Transfers

            //- Make changes following a parallel transfer
            template<class TrackCloudType>
            void correctAfterParallelTransfer(TrackCloudType&, trackingData&);


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            template<class TrackCloudType>
            bool hitPatch(TrackCloudType& cloud, trackingData& td);

            //- Overridable function to handle the particle hitting a wedgePatch
            template<class TrackCloudType>
            void hitWedgePatch(TrackCloudType&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  symmetryPlanePatch
            template<class TrackCloudType>
            void hitSymmetryPlanePatch(TrackCloudType&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  symmetryPatch
            template<class TrackCloudType>
            void hitSymmetryPatch(TrackCloudType&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  cyclicPatch
            template<class TrackCloudType>
            void hitCyclicPatch(TrackCloudType&, trackingData&);

            //- Overridable function to handle the particle hitting an
            //  nonConformalCyclicPolyPatch
            template<class TrackCloudType>
            bool hitNonConformalCyclicPatch
            (
                const vector& displacement,
                const scalar fraction,
                const label patchi,
                TrackCloudType& cloud,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            using ParcelType::hitProcessorPatch;

            //- Overridable function to handle the particle hitting a wallPatch
            template<class TrackCloudType>
            void hitWallPatch(TrackCloudType& cloud, trackingData& td);

            //- Overridable function to handle the particle hitting a basic
            //  patch. Fall-through for the above.
            template<class TrackCloudType>
            void hitBasicPatch(TrackCloudType& cloud, trackingData& td);


        // I-O

            //- Read
            template<class TrackCloudType>
            static void readFields(TrackCloudType& c);

            //- Write
            template<class TrackCloudType>
            static void writeFields(const TrackCloudType& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const MomentumParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MomentumParcelI.H"
#include "MomentumParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MomentumParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
