/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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
    Foam::ParticleCollector

Description
    Function object to collect the parcel mass- and mass flow rate over a
    set of polygons.  The polygons can either be specified by sets of user-
    supplied points, or in a concentric circles arrangement.  If a
    parcel is 'collected', it can be flagged to be removed from the
    domain using the removeCollected entry.

    Example usage:
    \verbatim
    particleCollector1
    {
        type            particleCollector;

        mode            concentricCircle;
        origin          (0.05 0.025 0.005);
        radius          (0.01 0.025 0.05);
        nSector         10;
        refDir          (1 0 0);
        normal          (0 0 1);

        negateParcelsOppositeNormal no;
        removeCollected no;
        surfaceFormat   vtk;
        resetOnWrite    no;
        log             yes;
    }

    particleCollector2
    {
        type            particleCollector;

        mode            polygon;
        polygons
        (
            (
                (0 0 0)
                (1 0 0)
                (1 1 0)
                (0 1 0)
            )
            (
                (0 0 1)
                (1 0 1)
                (1 1 1)
                (0 1 1)
            )
        );
        normal          (0 0 1);

        negateParcelsOppositeNormal no;
        removeCollected no;
        surfaceFormat   vtk;
        resetOnWrite    no;
        log             yes;
    }
    \endverbatim

SourceFiles
    ParticleCollector.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleCollector_H
#define ParticleCollector_H

#include "CloudFunctionObject.H"
#include "cylindricalCS.H"
#include "face.H"
#include "Switch.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ParticleCollector Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleCollector
:
    public CloudFunctionObject<CloudType>
{
public:

    enum modeType
    {
        mtPolygon,
        mtConcentricCircle,
        mtUnknown
    };


private:

    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;

        //- Collector mode type
        modeType mode_;

        //- Index of parcel types to collect (-1 by default = all particles)
        const label parcelType_;

        //- Flag to remove collected particles
        Switch removeCollected_;

        //- List of points
        Field<point> points_;

        //- List of faces
        List<face> faces_;


        // Polygon collector

            //- Triangulation of faces
            List<List<face> > faceTris_;

        // Concentric circles collector

            //- Number of sectors per circle
            label nSector_;

            //- List of radii
            List<scalar> radius_;

            //- Cylindrical co-ordinate system
            cylindricalCS coordSys_;


        //- Face areas
        Field<scalar> area_;

        //- Polygon normal vector per face
        Field<vector> normal_;

        //- Remove mass of parcel travelling in opposite direction to normal_
        bool negateParcelsOppositeNormal_;

        //- Surface output format
        const word surfaceFormat_;

        //- Flag to indicate whether data should be reset/cleared on writing
        Switch resetOnWrite_;

        //- Total time
        scalar totalTime_;

        //- Mass storage
        List<scalar> mass_;

        //- Mass total storage
        List<scalar> massTotal_;

        //- Mass flow rate storage
        List<scalar> massFlowRate_;

        //- Flag to indicate whether data should be written to file
        Switch log_;

        //- Output file pointer
        autoPtr<OFstream> outputFilePtr_;

        //- Last calculation time
        scalar timeOld_;

        //- Work list to store which faces are hit
        mutable DynamicList<label> hitFaceIDs_;


    // Private Member Functions

        //- Helper function to create log files
        void makeLogFile
        (
            const faceList& faces,
            const Field<point>& points,
            const Field<scalar>& area
        );

        //- Initialise polygon collectors
        void initPolygons(const List<Field<point> >& polygons);

        //- Initialise concentric circle collectors
        void initConcentricCircles();

        //- Collect parcels in polygon collectors
        void collectParcelPolygon
        (
            const point& p1,
            const point& p2
        ) const;

        //- Collect parcels in concentric circle collectors
        void collectParcelConcentricCircles
        (
            const point& p1,
            const point& p2
        ) const;


protected:

    // Protected Member Functions

        //- Write post-processing info
        void write();


public:

    //- Runtime type information
    TypeName("particleCollector");


    // Constructors

        //- Construct from dictionary
        ParticleCollector
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ParticleCollector(const ParticleCollector<CloudType>& pc);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType> > clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType> >
            (
                new ParticleCollector<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleCollector();


    // Member Functions

        // Access

            //- Return const access to the reset on write flag
            inline const Switch& resetOnWrite() const;


        // Evaluation

            //- Post-move hook
            virtual void postMove
            (
                parcelType& p,
                const label cellI,
                const scalar dt,
                const point& position0,
                bool& keepParticle
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleCollectorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ParticleCollector.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
