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
    Foam::star

Description
    Engine for constructing a star-shaped domain by walking

SourceFiles
    star.C

\*---------------------------------------------------------------------------*/

#ifndef star_H
#define star_H

#include "DynamicList.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class star Declaration
\*---------------------------------------------------------------------------*/

class star
{
public:

    // Public Classes

        //- Star-edge structure. Poor man's linked list link.
        struct edgeLink
        {
            label starEdgei0_;
            label edgei_;
            label starEdgei1_;
        };

        //- Context class. Returned by populate and resets the star when it
        //  goes out of scope.
        class context
        {
        private:

            // Private Data

                //- Reference to the star
                star& star_;


        public:

            // Constructors

                //- Construct from components
                context(star& s);


            // Destructor
            ~context();
        };


private:

    // Private Data

        //- For each star face, the external face index, or -1
        DynamicList<label> starFaceFaces_;

        //- For each external face, the star face index, or -1
        DynamicList<label> faceStarFaces_;

        //- For each star edge, the external edge index and the two adjacent
        //  star indices, or {-1, -1, -1}
        DynamicList<edgeLink> starEdgeEdges_;

        //- For each external edge, the star edge index, or -1
        DynamicList<label> edgeStarEdges_;

        //- Generic work array
        DynamicList<label> work_;


    // Private Member Functions

        //- Reset the contents of the star, ready for another populate
        void reset();

        //- Swap the given star edges. The star connectivity remains the same.
        //  Used to change order of iteration.
        void swapStarEdges(const label starEdgeiA, const label starEdgeiB);

        //- Expand the star by crossing a given edge into a given face
        template<class FaceEdges>
        void cross
        (
            const label edgei,
            const label facei,
            const UList<FaceEdges>& faceEdges
        );


public:

    // Constructors

        //- Construct null
        star();


    // Destructor
    ~star();


    // Member Functions

        // Access

            //- Access the star-face-faces
            const DynamicList<label>& starFaceFaces() const
            {
                return starFaceFaces_;
            }

            //- Access the face-star-faces
            const DynamicList<label>& faceStarFaces() const
            {
                return faceStarFaces_;
            }

            //- Access the star-edge-edges
            const DynamicList<edgeLink>& starEdgeEdges() const
            {
                return starEdgeEdges_;
            }

            //- Access the edge-star-edges
            const DynamicList<label>& edgeStarEdges() const
            {
                return edgeStarEdges_;
            }


        //- Populate the star given a seed face or edge, a function
        //  `canCross(edgei, facei)` to determine whether an edge may be cross,
        //  and face-edge and edge-face addressing. After this is called, the
        //  forAll... macros below can be user to iterate through the star
        //  faces, or around the edges of the star perimeter.
        template<class CanCross, class FaceEdges>
        context populate
        (
            const label faceOrEdgei,
            const bool isFace,
            CanCross canCross,
            const UList<FaceEdges>& faceEdges,
            const UList<labelPair>& edgeFaces,
            const label maxIterations = labelMax
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline bool operator==(const star::edgeLink& elA, const star::edgeLink& elB)
{
    return
        elA.starEdgei0_ == elB.starEdgei0_
     && elA.edgei_ == elB.edgei_
     && elA.starEdgei1_ == elB.starEdgei1_;
}


inline bool operator!=(const star::edgeLink& elA, const star::edgeLink& elB)
{
    return !(elA == elB);
}


inline Ostream& operator<<(Ostream& os, const star::edgeLink& el)
{
    return os
        << '(' << el.starEdgei0_ << ')'
        << ' ' << el.edgei_ << ' '
        << '(' << el.starEdgei1_ << ')';
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forAllStarFaces(star, starFacei, facei)                                \
    for                                                                        \
    (                                                                          \
        label starFacei = 0,                                                   \
        facei =                                                                \
        (                                                                      \
            star.starFaceFaces().size() > 0                                    \
          ? (star).starFaceFaces()[starFacei]                                  \
          : -1                                                                 \
        );                                                                     \
        starFacei < star.starFaceFaces().size();                               \
        ++ starFacei,                                                          \
        facei =                                                                \
        (                                                                      \
            starFacei < (star).starFaceFaces().size()                          \
          ? (star).starFaceFaces()[starFacei]                                  \
          : -1                                                                 \
        )                                                                      \
    )


#define forAllStarEdges(star, i, starEdgei, edgei)                             \
    for                                                                        \
    (                                                                          \
        label i = 0,                                                           \
        starEdgei = 0,                                                         \
        edgei =                                                                \
        (                                                                      \
            (star).starEdgeEdges().size() > 0                                  \
          ? (star).starEdgeEdges()[0].edgei_                                   \
          : -1                                                                 \
        );                                                                     \
        i < (star).starEdgeEdges().size();                                     \
        ++ i,                                                                  \
        starEdgei = (star).starEdgeEdges()[starEdgei].starEdgei1_,             \
        edgei =                                                                \
        (                                                                      \
            starEdgei != -1                                                    \
          ? (star).starEdgeEdges()[starEdgei].edgei_                           \
          : -1                                                                 \
        )                                                                      \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "starTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
