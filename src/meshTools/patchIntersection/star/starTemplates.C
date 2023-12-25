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

\*---------------------------------------------------------------------------*/

#include "OFstream.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceEdges>
void Foam::star::cross
(
    const label edgei,
    const label facei,
    const UList<FaceEdges>& faceEdges
)
{
    // Add the new face to the star
    starFaceFaces_.append(facei);
    faceStarFaces_[facei] = starFaceFaces_.size() - 1;

    // Store information about this edge, then remove it from the star
    const label faceEdgei = findIndex(faceEdges[facei], edgei);
    const label starEdgei = edgeStarEdges_[faceEdges[facei][faceEdgei]];
    const label starEdgei0 = starEdgeEdges_[starEdgei].starEdgei0_;
    const label starEdgei1 = starEdgeEdges_[starEdgei].starEdgei1_;
    edgeStarEdges_[starEdgeEdges_[starEdgei].edgei_] = -1;
    starEdgeEdges_[starEdgei] = {-1, -1, -1};

    // Walk backwards and forwards around the face until we find edges not
    // already in the star. Remove edges in the star as we go.
    label faceEdgeiR = faceEdgei, starEdgeiR1 = starEdgei1;
    while (true)
    {
        const label faceEdgeiRR = faceEdges[facei].rcIndex(faceEdgeiR);
        const label edgeiRR = faceEdges[facei][faceEdgeiRR];
        const label starEdgeiRR = edgeStarEdges_[edgeiRR];
        if (starEdgeiRR == -1 || starEdgeiRR != starEdgeiR1) break;
        faceEdgeiR = faceEdgeiRR;
        starEdgeiR1 = starEdgeEdges_[starEdgeiRR].starEdgei1_;
        edgeStarEdges_[starEdgeEdges_[starEdgeiRR].edgei_] = -1;
        starEdgeEdges_[starEdgeiRR] = {-1, -1, -1};
    }
    label faceEdgeiF = faceEdgei, starEdgeiF0 = starEdgei0;
    while (true)
    {
        const label faceEdgeiFF = faceEdges[facei].fcIndex(faceEdgeiF);
        const label edgeiFF = faceEdges[facei][faceEdgeiFF];
        const label starEdgeiFF = edgeStarEdges_[edgeiFF];
        if (starEdgeiFF == -1 || starEdgeiFF != starEdgeiF0) break;
        faceEdgeiF = faceEdgeiFF;
        starEdgeiF0 = starEdgeEdges_[starEdgeiFF].starEdgei0_;
        edgeStarEdges_[starEdgeEdges_[starEdgeiFF].edgei_] = -1;
        starEdgeEdges_[starEdgeiFF] = {-1, -1, -1};
    }

    // Get the first edge after the forwards edge
    label faceEdgej = faceEdges[facei].fcIndex(faceEdgeiF);

    // If there are no face edges not yet in the star, then all edges are
    // to be removed. Just connect up the adjacent edges.
    if (faceEdgej == faceEdgeiR)
    {
        starEdgeEdges_[starEdgeiF0].starEdgei1_ = starEdgeiR1;
        starEdgeEdges_[starEdgeiR1].starEdgei0_ = starEdgeiF0;
    }

    // If there are face edges not yet in the star, then loop over them and
    // add them into the star
    while (faceEdgej != faceEdgeiR)
    {
        const bool isFirst =
            faceEdgej == faceEdges[facei].fcIndex(faceEdgeiF);
        const bool isLast =
            faceEdgej == faceEdges[facei].rcIndex(faceEdgeiR);

        const label starEdgej0 =
            isFirst ? starEdgeiF0 : starEdgeEdges_.size() - 1;
        const label starEdgej = starEdgeEdges_.size();
        const label starEdgej1 =
            isLast ? starEdgeiR1 : starEdgeEdges_.size() + 1;

        const label edgej = faceEdges[facei][faceEdgej];

        if (isFirst && starEdgeiF0 != -1)
        {
            starEdgeEdges_[starEdgeiF0].starEdgei1_ = starEdgej;
        }

        starEdgeEdges_.append({starEdgej0, edgej, starEdgej1});
        edgeStarEdges_[edgej] = starEdgej;

        if (isLast && starEdgeiR1 != -1)
        {
            starEdgeEdges_[starEdgeiR1].starEdgei0_ = starEdgej;
        }

        faceEdgej = faceEdges[facei].fcIndex(faceEdgej);
    }
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CanCross, class FaceEdges>
Foam::star::context Foam::star::populate
(
    const label faceOrEdgei,
    const bool isFace,
    CanCross canCross,
    const UList<FaceEdges>& faceEdges,
    const UList<labelPair>& edgeFaces,
    const label maxIterations
)
{
    // Expand to fit supplied data
    faceStarFaces_.resize(faceEdges.size(), -1);
    edgeStarEdges_.resize(edgeFaces.size(), -1);

    // Initialise
    if (isFace)
    {
        // Add one edge from the to the star and then cross into the face
        const label edgei = faceEdges[faceOrEdgei][0];
        starEdgeEdges_.append({-1, edgei, -1});
        edgeStarEdges_[edgei] = 0;
        cross(edgei, faceOrEdgei, faceEdges);

        // Add the edge back again as it was removed by the crossing. Reuse the
        // first element of the starEdgeEdges array. Join up to close the star.
        const label starEdgei0 = starEdgeEdges_.size() - 1;
        const label starEdgei1 = 1;
        starEdgeEdges_[0] = {starEdgei0, edgei, starEdgei1};
        edgeStarEdges_[edgei] = 0;
        starEdgeEdges_[starEdgei0].starEdgei1_ = 0;
        starEdgeEdges_[starEdgei1].starEdgei0_ = 0;
    }
    else
    {
        bool first = true;

        // Loop the adjacent faces
        forAll(edgeFaces[faceOrEdgei], edgeFacei)
        {
            const label facei = edgeFaces[faceOrEdgei][edgeFacei];
            if (facei == -1) continue;

            // If the first face, add the edge to the star
            if (first)
            {
                starEdgeEdges_.append({-1, faceOrEdgei, -1});
                edgeStarEdges_[faceOrEdgei] = 0;
                first = false;
            }
            // If the second face, add the edge into the star again as it was
            // removed by the previous crossing. Reuse the first element of the
            // starEdgeEdges array. Join up to close the star.
            else
            {
                const label starEdgei0 = starEdgeEdges_.size() - 1;
                const label starEdgei1 = 1;
                starEdgeEdges_[0] = {starEdgei0, faceOrEdgei, starEdgei1};
                edgeStarEdges_[faceOrEdgei] = 0;
                starEdgeEdges_[starEdgei0].starEdgei1_ = 0;
                starEdgeEdges_[starEdgei1].starEdgei0_ = 0;
            }

            // Cross the edge into the face
            cross(faceOrEdgei, facei, faceEdges);
        }
    }

    // Walk the star to get intersected faces and their perimeter
    label iterations = 0;
    forAll(starEdgeEdges_, starEdgei)
    {
        if (iterations == maxIterations) break;
        ++ iterations;

        const label edgei = starEdgeEdges_[starEdgei].edgei_;
        if (edgei == -1) continue;

        // Get the adjacent non-star face
        label facei = -1;
        forAll(edgeFaces[edgei], edgeFacei)
        {
            facei = edgeFaces[edgei][edgeFacei];
            if (facei != -1 && faceStarFaces_[facei] == -1)
            {
                break;
            }
            facei = -1;
        }

        // If the crossing condition is satisfied then expand the star
        if (facei != -1 && canCross(edgei, facei))
        {
            cross(edgei, facei, faceEdges);
        }
    }

    // Remove spikes
    forAll(starEdgeEdges_, starEdgei)
    {
        label starEdgei0 = starEdgei;
        label starEdgei1 = starEdgeEdges_[starEdgei0].starEdgei1_;

        if (starEdgei1 == -1) continue;

        label edgei0 = starEdgeEdges_[starEdgei0].edgei_;
        label edgei1 = starEdgeEdges_[starEdgei1].edgei_;

        while (edgei0 != -1 && edgei1 != -1 && edgei0 == edgei1)
        {
            const label starEdgei00 = starEdgeEdges_[starEdgei0].starEdgei0_;
            const label starEdgei11 = starEdgeEdges_[starEdgei1].starEdgei1_;

            edgeStarEdges_[edgei0] = -1;
            starEdgeEdges_[starEdgei0] = {-1, -1, -1};
            edgeStarEdges_[edgei1] = -1;
            starEdgeEdges_[starEdgei1] = {-1, -1, -1};

            if (starEdgei00 != -1)
            {
                starEdgeEdges_[starEdgei00].starEdgei1_ = starEdgei11;
                starEdgei0 = starEdgei00;
                edgei0 = starEdgeEdges_[starEdgei00].edgei_;
            }
            else
            {
                starEdgei0 = edgei0 = -1;
            }

            if (starEdgei11 != -1)
            {
                starEdgeEdges_[starEdgei11].starEdgei0_ = starEdgei00;
                starEdgei1 = starEdgei11;
                edgei1 = starEdgeEdges_[starEdgei11].edgei_;
            }
            else
            {
                starEdgei1 = edgei0 = -1;
            }
        }
    }

    // Allocate map storage
    DynamicList<label>& oldStarEdgeNewStarEdges = work_;
    oldStarEdgeNewStarEdges.resize(starEdgeEdges_.size());
    oldStarEdgeNewStarEdges = -1;

    // Remove empty edges by shuffling up
    label starEdgei = 0;
    forAll(starEdgeEdges_, starEdgej)
    {
        if (starEdgeEdges_[starEdgej].edgei_ != -1)
        {
            oldStarEdgeNewStarEdges[starEdgej] = starEdgei;
            starEdgeEdges_[starEdgei] = starEdgeEdges_[starEdgej];
            ++ starEdgei;
        }
    }
    starEdgeEdges_.resize(starEdgei);

    // Map
    forAll(starEdgeEdges_, starEdgei)
    {
        label& starEdgei0 = starEdgeEdges_[starEdgei].starEdgei0_;
        if (starEdgei0 != -1)
        {
            starEdgei0 = oldStarEdgeNewStarEdges[starEdgei0];
        }
        label& starEdgei1 = starEdgeEdges_[starEdgei].starEdgei1_;
        if (starEdgei1 != -1)
        {
            starEdgei1 = oldStarEdgeNewStarEdges[starEdgei1];
        }
    }

    /*
    // Randomise swapping
    static Random rndGen(12345);
    const label starEdgeiA = rndGen.sampleAB<label>(0, starEdgeEdges_.size());
    const label starEdgeiB = rndGen.sampleAB<label>(0, starEdgeEdges_.size());
    swapStarEdges(starEdgeiA, starEdgeiB);
    */

    // If the star is open, put the first edge first, so that forAllStarEdges
    // always starts in the right place
    forAll(starEdgeEdges_, starEdgei)
    {
        if (starEdgeEdges_[starEdgei].starEdgei0_ == -1)
        {
            swapStarEdges(starEdgei, 0);
        }
    }

    return context(*this);
}


// ************************************************************************* //
