/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::KinematicParcel

Description
    Kinematic parcel class with rotational motion (as spherical
    particles only) and one/two-way coupling with the continuous
    phase.

    Sub-models include:
    - drag
    - turbulent dispersion
    - wall interactions

SourceFiles
    KinematicParcelI.H
    KinematicParcel.C
    KinematicParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef KinematicParcel_H
#define KinematicParcel_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "interpolation.H"
#include "demandDrivenEntry.H"

// #include "ParticleForceList.H" // TODO

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class KinematicParcel;

// Forward declaration of friend functions

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const KinematicParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                         Class KinematicParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class KinematicParcel
:
    public ParcelType
{
    // Private data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;

        //- Number of particle tracking attempts before we assume that it stalls
        static label maxTrackAttempts;


public:

    //- Class to hold kinematic particle constant properties
    class constantProperties
    {
    protected:

        // Protected data

            //- Constant properties dictionary
            const dictionary dict_;


    private:

        // Private data

            //- Parcel type id - used for post-processing to flag the type
            //  of parcels issued by this cloud
            demandDrivenEntry<label> parcelTypeId_;

            //- Minimum density [kg/m3]
            demandDrivenEntry<scalar> rhoMin_;

            //- Particle density [kg/m3] (constant)
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


        // Member functions

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
    public:

        enum trackPart
        {
            tpVelocityHalfStep,
            tpLinearTrack,
            tpRotationalTrack
        };


    private:

        // Private data

            // Interpolators for continuous phase fields

                //- Density interpolator
                autoPtr<interpolation<scalar>> rhoInterp_;

                //- Velocity interpolator
                autoPtr<interpolation<vector>> UInterp_;

                //- Dynamic viscosity interpolator
                autoPtr<interpolation<scalar>> muInterp_;


            // Cached continuous phase properties

                //- Density [kg/m3]
                scalar rhoc_;

                //- Velocity [m/s]
                vector Uc_;

                //- Viscosity [Pa.s]
                scalar muc_;


            //- Local gravitational or other body-force acceleration
            const vector& g_;

            // label specifying which part of the integration
            // algorithm is taking place
            trackPart part_;


    public:

        // Constructors

            //- Construct from components
            template <class TrackCloudType>
            inline trackingData
            (
                const TrackCloudType& cloud,
                trackPart part = tpLinearTrack
            );


        // Member functions

            //- Return conat access to the interpolator for continuous
            //  phase density field
            inline const interpolation<scalar>& rhoInterp() const;

            //- Return conat access to the interpolator for continuous
            //  phase velocity field
            inline const interpolation<vector>& UInterp() const;

            //- Return conat access to the interpolator for continuous
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

            // Return const access to the gravitational acceleration vector
            inline const vector& g() const;

            //- Return the part of the tracking operation taking place
            inline trackPart part() const;

            //- Return access to the part of the tracking operation taking place
            inline trackPart& part();
    };


protected:

    // Protected data

        // Parcel properties

            //- Active flag - tracking inactive when active = false
            bool active_;

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

            //- Density [kg/m3]
            scalar rho_;

            //- Age [s]
            scalar age_;

            //- Time spent in turbulent eddy [s]
            scalar tTurb_;

            //- Turbulent velocity fluctuation [m/s]
            vector UTurb_;


    // Protected Member Functions

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

    // Static data members

        //- Runtime type information
        TypeName("KinematicParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " active"
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
        inline KinematicParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline KinematicParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline KinematicParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const scalar dTarget0,
            const vector& U0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        KinematicParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        KinematicParcel(const KinematicParcel& p);

        //- Construct as a copy
        KinematicParcel(const KinematicParcel& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new KinematicParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new KinematicParcel(*this, mesh));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<KinematicParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<KinematicParcel<ParcelType>>
                (
                    new KinematicParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to active flag
            inline bool active() const;

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

            //- Return const access to active flag
            inline bool& active();

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


        // Main calculation loop

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


        // Tracking

            //- Move the parcel
            template<class TrackCloudType>
            bool move
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar trackTime
            );


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            template<class TrackCloudType>
            bool hitPatch(TrackCloudType& cloud, trackingData& td);

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            template<class TrackCloudType>
            void hitProcessorPatch(TrackCloudType& cloud, trackingData& td);

            //- Overridable function to handle the particle hitting a wallPatch
            template<class TrackCloudType>
            void hitWallPatch(TrackCloudType& cloud, trackingData& td);

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            virtual void transformProperties(const vector& separation);


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
            const KinematicParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "KinematicParcelI.H"
#include "KinematicParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "KinematicParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
