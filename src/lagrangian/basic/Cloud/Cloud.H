/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::Cloud

Description
    Base cloud calls templated on particle type

SourceFiles
    Cloud.C
    CloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef Cloud_H
#define Cloud_H

#include "cloud.H"
#include "IDLList.H"
#include "IOField.H"
#include "CompactIOField.H"
#include "polyMesh.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of functions
template<class ParticleType>
class Cloud;

template<class ParticleType>
class IOPosition;

template<class ParticleType>
Ostream& operator<<
(
    Ostream&,
    const Cloud<ParticleType>&
);


/*---------------------------------------------------------------------------*\
                           Class Cloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class Cloud
:
    public cloud,
    public IDLList<ParticleType>
{
    // Private data

        const polyMesh& polyMesh_;

        //- Temporary storage for addressing. Used in findTris.
        mutable DynamicList<label> labels_;

        //- Count of how many tracking rescue corrections have been
        //  applied
        mutable label nTrackingRescues_;

        //- Does the cell have wall faces
        mutable autoPtr<PackedBoolList> cellWallFacesPtr_;


    // Private Member Functions

        //- Check patches
        void checkPatches() const;

        //- Initialise cloud on IO constructor
        void initCloud(const bool checkClass);

        //- Find all cells which have wall faces
        void calcCellWallFaces() const;

        //- Read cloud properties dictionary
        void readCloudUniformProperties();

        //- Write cloud properties dictionary
        void writeCloudUniformProperties() const;


public:

    friend class particle;
    template<class ParticleT>
    friend class IOPosition;

    typedef ParticleType particleType;

    typedef typename IDLList<ParticleType>::iterator iterator;
    typedef typename IDLList<ParticleType>::const_iterator const_iterator;

    //-Runtime type information
    TypeName("Cloud");


    // Static data

        //- Name of cloud properties dictionary
        static word cloudPropertiesName;


    // Constructors

        //- Construct from mesh and a list of particles
        Cloud
        (
            const polyMesh& mesh,
            const IDLList<ParticleType>& particles
        );

        //- Construct from mesh, cloud name, and a list of particles
        Cloud
        (
            const polyMesh& mesh,
            const word& cloudName,
            const IDLList<ParticleType>& particles
        );

        //- Construct from mesh by reading from file
        //  Optionally disable checking of class name for post-processing
        Cloud
        (
            const polyMesh& mesh,
            const bool checkClass = true
        );


        //- Construct from mesh by reading from file with given cloud instance
        //  Optionally disable checking of class name for post-processing
        Cloud
        (
            const polyMesh& pMesh,
            const word& cloudName,
            const bool checkClass = true
        );


    // Member Functions

        // Access

            //- Return the polyMesh reference
            const polyMesh& pMesh() const
            {
                return polyMesh_;
            }

            label size() const
            {
                return IDLList<ParticleType>::size();
            };

            DynamicList<label>& labels()
            {
                return labels_;
            }

            //- Return nTrackingRescues
            label nTrackingRescues() const
            {
                return nTrackingRescues_;
            }

            //- Increment the nTrackingRescues counter
            void trackingRescue() const
            {
                nTrackingRescues_++;
                if (cloud::debug && size() && (nTrackingRescues_ % size() == 0))
                {
                    Pout<< "    " << nTrackingRescues_
                        << " tracking rescues " << endl;
                }
            }

            //- Whether each cell has any wall faces (demand driven data)
            const PackedBoolList& cellHasWallFaces() const;

            //- Switch to specify if particles of the cloud can return
            //  non-zero wall distance values.  By default, assume
            //  that they can't (default for wallImpactDistance in
            //  particle is 0.0).
            virtual bool hasWallImpactDistance() const
            {
                return false;
            }


            // Iterators

                const const_iterator begin() const
                {
                    return IDLList<ParticleType>::begin();
                };

                const const_iterator cbegin() const
                {
                    return IDLList<ParticleType>::cbegin();
                };

                const const_iterator end() const
                {
                    return IDLList<ParticleType>::end();
                };

                const const_iterator cend() const
                {
                    return IDLList<ParticleType>::cend();
                };

                iterator begin()
                {
                    return IDLList<ParticleType>::begin();
                };

                iterator end()
                {
                    return IDLList<ParticleType>::end();
                };


        // Edit

            void clear()
            {
                IDLList<ParticleType>::clear();
            };

            //- Transfer particle to cloud
            void addParticle(ParticleType* pPtr);

            //- Remove particle from cloud and delete
            void deleteParticle(ParticleType&);

            //- Reset the particles
            void cloudReset(const Cloud<ParticleType>& c);

            //- Move the particles
            //  passing the TrackingData to the track function
            template<class TrackData>
            void move(TrackData& td, const scalar trackTime);

            //- Remap the cells of particles corresponding to the
            //  mesh topology change
            template<class TrackData>
            void autoMap(TrackData& td, const mapPolyMesh&);


        // Read

            //- Helper to construct IOobject for field and current time.
            IOobject fieldIOobject
            (
                const word& fieldName,
                const IOobject::readOption r
            ) const;

            //- Check lagrangian data field
            template<class DataType>
            void checkFieldIOobject
            (
                const Cloud<ParticleType>& c,
                const IOField<DataType>& data
            ) const;

            //- Check lagrangian data fieldfield
            template<class DataType>
            void checkFieldFieldIOobject
            (
                const Cloud<ParticleType>& c,
                const CompactIOField<Field<DataType>, DataType>& data
            ) const;

            //- Read the field data for the cloud of particles. Dummy at
            //  this level.
            virtual void readFields();


        // Write

            //- Write the field data for the cloud of particles Dummy at
            //  this level.
            virtual void writeFields() const;

            //- Write using given format, version and compression.
            //  Only writes the cloud file if the Cloud isn't empty
            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp
            ) const;

            //- Write positions to \<cloudName\>_positions.obj file
            void writePositions() const;


    // Ostream Operator

        friend Ostream& operator<< <ParticleType>
        (
            Ostream&,
            const Cloud<ParticleType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Cloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
