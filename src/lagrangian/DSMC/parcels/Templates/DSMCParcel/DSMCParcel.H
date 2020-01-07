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
    Foam::DSMCParcel

Description
    DSMC parcel class

SourceFiles
    DSMCParcelI.H
    DSMCParcel.C
    DSMCParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef DSMCParcel_H
#define DSMCParcel_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "contiguous.H"
#include "DSMCCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class DSMCParcel;

// Forward declaration of friend functions

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const DSMCParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                         Class DSMCParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class DSMCParcel
:
    public ParcelType
{
    // Private member data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


public:

    //- Class to hold DSMC particle constant properties
    class constantProperties
    {
        // Private Data

            //- Particle mass [kg] (constant)
            scalar mass_;

            //- Particle hard sphere diameter [m] (constant)
            scalar d_;

            //- Internal degrees of freedom
            direction internalDegreesOfFreedom_;

            //- Viscosity index
            scalar omega_;


    public:

        // Constructors

            //- Null constructor, allows List of constantProperties to be
            //  created before the contents is initialised
            inline constantProperties();

            //- Constructor from dictionary
            inline constantProperties(const dictionary& dict);


        // Member Functions

            //- Return const access to the particle mass [kg]
            inline scalar mass() const;

            //- Return const access to the hard sphere diameter [m]
            inline scalar d() const;

            //- Return the reference total collision cross section
            inline scalar sigmaT() const;

            //- Return the internalDegreesOfFreedom
            inline direction internalDegreesOfFreedom() const;

            //- Return the viscosity index
            inline scalar omega() const;

    };


    //- Use base tracking data
    typedef typename ParcelType::trackingData trackingData;


protected:

    // Protected member data

        // Parcel properties

            //- Velocity of Parcel [m/s]
            vector U_;

            //- Internal energy of the Parcel, covering all non-translational
            //  degrees of freedom [J]
            scalar Ei_;

            //- Parcel type id
            label typeId_;


public:

    //- Runtime type information
    TypeName("DSMCParcel");

    friend class Cloud<ParcelType>;


    // Constructors

        //- Construct from components
        inline DSMCParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const vector& U,
            const scalar Ei,
            const label typeId
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology
        inline DSMCParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const vector& U,
            const scalar Ei,
            const label typeId
        );

        //- Construct from Istream
        DSMCParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new DSMCParcel<ParcelType>(*this));
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

            autoPtr<DSMCParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<DSMCParcel<ParcelType>>
                (
                    new DSMCParcel<ParcelType>(mesh_, is, true)
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
            bool hitPatch(TrackCloudType&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            template<class TrackCloudType>
            void hitProcessorPatch(TrackCloudType&, trackingData&);

            //- Overridable function to handle the particle hitting a wallPatch
            template<class TrackCloudType>
            void hitWallPatch(TrackCloudType&, trackingData&);

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const transformer&);


        // I-O

            static void readFields(Cloud<DSMCParcel<ParcelType>>& c);

            static void writeFields(const Cloud<DSMCParcel<ParcelType>>& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const DSMCParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DSMCParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DSMCParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
