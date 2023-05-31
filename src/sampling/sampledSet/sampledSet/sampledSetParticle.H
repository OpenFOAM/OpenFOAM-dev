/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    Foam::sampledSetParticle

Description
    Particle for generating line-type sampled sets

SourceFiles
    sampledSetParticle.C

\*---------------------------------------------------------------------------*/

#ifndef sampledSetParticle_H
#define sampledSetParticle_H

#include "particle.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class sampledSetParticle;
class sampledSetCloud;

Ostream& operator<<(Ostream&, const sampledSetParticle&);


/*---------------------------------------------------------------------------*\
                          Class sampledSetParticle Declaration
\*---------------------------------------------------------------------------*/

class sampledSetParticle
:
    public particle
{
public:

    //- Tracking data class
    class trackingData
    :
        public particle::trackingData
    {
    public:

        //- List of set points to track through
        const List<point>& set_;

        //- Are we storing the set points?
        const bool storeSet_;

        //- Are we storing the faces we track through?
        const label storeFaces_;

        //- Are we storing the cells we track through?
        const bool storeCells_;

        //- Stored positions
        DynamicList<point>& positions_;

        //- Stored distances
        DynamicList<scalar>& distances_;

        //- Stored cell indices
        DynamicList<label>& cells_;

        //- Stored face indices
        DynamicList<label>& faces_;

        //- Construct from components
        trackingData
        (
            sampledSetCloud& cloud,
            const List<point>& set,
            const bool storeSet,
            const label storeFaces,
            const bool storeCells,
            DynamicList<point>& positions,
            DynamicList<scalar>& distances,
            DynamicList<label>& cells,
            DynamicList<label>& faces
        )
        :
            particle::trackingData(cloud),
            set_(set),
            storeSet_(storeSet),
            storeFaces_(storeFaces),
            storeCells_(storeCells),
            positions_(positions),
            distances_(distances),
            cells_(cells),
            faces_(faces)
        {}
    };


private:

    // Private Data

        //- The set interval we are currently within
        label seti_;

        //- The fraction of this set interval from the end
        scalar setF_;

        //- The scalar distance through the set
        scalar distance_;

        //- Do we have the previous face position stored?
        bool havePosition0_;

        //- The previous face position
        point position0_;


public:

    // Constructors

        //- Construct from components
        sampledSetParticle
        (
            const polyMesh& mesh,
            const point& position,
            const label celli,
            const label seti,
            const scalar setF,
            const scalar distance
        );

        //- Construct from Istream
        sampledSetParticle(Istream& is, bool readFields = true);

        //- Copy constructor
        sampledSetParticle(const sampledSetParticle&);

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new sampledSetParticle(*this));
        }

        //- Construct from Istream and return
        static autoPtr<sampledSetParticle> New(Istream& is)
        {
            return autoPtr<sampledSetParticle>(new sampledSetParticle(is));
        }


    // Member Functions

        // Storage

            //- Store a point
            void store(sampledSetCloud&, trackingData&);

            //- Store a point on a face
            void storeFace(sampledSetCloud&, trackingData&);

            //- Store a point in a cell
            void storeCell(sampledSetCloud&, trackingData&);


        // Tracking

            //- Track all particles to their end point
            bool move(sampledSetCloud&, trackingData&);

            //- Hit a wedge patch. Ends the track.
            void hitWedgePatch(sampledSetCloud&, trackingData&);

            //- Hit a symmetry patch. Ends the track.
            void hitSymmetryPlanePatch(sampledSetCloud&, trackingData&);

            //- Hit a symmetry plane patch. Ends the track.
            void hitSymmetryPatch(sampledSetCloud&, trackingData&);

            //- Hit a cyclic patch. Ends the track.
            void hitCyclicPatch(sampledSetCloud&, trackingData&);

            //- Hit a processor patch. Transfers the particle if there is no
            //  transformation. Ends the track if there is a transformation.
            void hitProcessorPatch(sampledSetCloud&, trackingData&);

            //- Hit a wall patch. Ends the track.
            void hitWallPatch(sampledSetCloud&, trackingData&);

            //- Do corrections to the particle and tracking data following a
            //  transfer between processors
            void correctAfterParallelTransfer(sampledSetCloud&, trackingData&);


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const sampledSetParticle&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
