/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::DsmcParcel

Description
    DSMC parcel class

SourceFiles
    DsmcParcelI.H
    DsmcParcel.C
    DsmcParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef DsmcParcel_H
#define DsmcParcel_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "contiguous.H"

#include "DsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class DsmcParcel;

// Forward declaration of friend functions

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const DsmcParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                         Class DsmcParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class DsmcParcel
:
    public ParcelType
{
public:

    //- Class to hold DSMC particle constant properties
    class constantProperties
    {
        // Private data

            //- Particle mass [kg] (constant)
            scalar mass_;

            //- Particle hard sphere diameter [m] (constant)
            scalar d_;

            //- Internal degrees of freedom
            scalar internalDegreesOfFreedom_;

            //- Viscosity index
            scalar omega_;


    public:

        // Constrcutors

            //- Null constructor, allows List of constantProperties to be
            //  created before the contents is initialised
            inline constantProperties();

            //- Constructor from dictionary
            inline constantProperties(const dictionary& dict);


        // Member functions

            //- Return const access to the particle density
            inline scalar mass() const;

            //- Return const access to the minimum particle mass
            inline scalar d() const;

            //- Return the reference total collision cross section
            inline scalar sigmaT() const;

            //- Return the internalDegreesOfFreedom
            inline scalar internalDegreesOfFreedom() const;

            //- Return the viscosity index
            inline scalar omega() const;

    };


    //- Class used to pass kinematic tracking data to the trackToFace function
    class trackingData
    :
        public particle::TrackingData<DsmcCloud<DsmcParcel<ParcelType> > >
    {
    public:

        // Constructors

            //- Construct from components
            trackingData(DsmcCloud<DsmcParcel<ParcelType> >& cloud)
            :
                particle::TrackingData<DsmcCloud<DsmcParcel<ParcelType> > >
                (
                    cloud
                )
            {}
    };


protected:

    // Protected member data

        // Parcel properties

            //- Velocity of Parcel [m/s]
            vector U_;

            //- Internal energy of the Parcel, covering all non-translational
            // degrees of freedom [J]
            scalar Ei_;

            //- Parcel type id
            label typeId_;


public:

    //- Runtime type information
    TypeName("DsmcParcel");

    friend class Cloud<ParcelType>;


    // Constructors

        //- Construct from components
        inline DsmcParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const vector& U,
            const scalar Ei,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const label typeId
        );

        //- Construct from Istream
        DsmcParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new DsmcParcel<ParcelType>(*this));
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

            autoPtr<DsmcParcel<ParcelType> > operator()(Istream& is) const
            {
                return autoPtr<DsmcParcel<ParcelType> >
                (
                    new DsmcParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return type id
            inline label typeId() const;

            //- Return const access to velocity
            inline const vector& U() const;

            //- Return const access to internal energy
            inline scalar Ei() const;

        // Edit

            //- Return access to velocity
            inline vector& U();

            //- Return access to internal energy
            inline scalar& Ei();


        // Main calculation loop

            // Tracking

                //- Move the parcel
                template<class TrackData>
                bool move(TrackData& td, const scalar trackTime);


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            template<class TrackData>
            bool hitPatch
            (
                const polyPatch&,
                TrackData& td,
                const label patchI,
                const scalar trackFraction,
                const tetIndices& tetIs
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            template<class TrackData>
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                TrackData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            template<class TrackData>
            void hitWallPatch
            (
                const wallPolyPatch&,
                TrackData& td,
                const tetIndices&
            );

            //- Overridable function to handle the particle hitting a polyPatch
            template<class TrackData>
            void hitPatch
            (
                const polyPatch&,
                TrackData& td
            );

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            virtual void transformProperties(const vector& separation);


        // I-O

            static void readFields(Cloud<DsmcParcel<ParcelType> >& c);

            static void writeFields(const Cloud<DsmcParcel<ParcelType> >& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const DsmcParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DsmcParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DsmcParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
