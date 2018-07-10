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
    Foam::topoCellLooper

Description
    Implementation of cellLooper. This one recognizes splitHexes and tries
    to make a cut such that if the neighbour was split (in a previous iteration)
    this one also gets split in the same direction so that the result
    will be a mesh without splitHexes.

    'splitHexes' are cells of which the 'featureEdges'
    (see cellFeatures class) form a hex. The remaining non-feature edges
    are assumed to result from splitting the neighbour and this class tries
    to start from one of these and cut through to an opposite edge.

    The current set of cuts (vertIsCut, edgeIsCut, edgeWeight) are not being
    used by this implementation.

    All non-splitHexes are done by the parent classes.


SourceFiles
    topoCellLooper.C

\*---------------------------------------------------------------------------*/

#ifndef topoCellLooper_H
#define topoCellLooper_H

#include "hexCellLooper.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class cellFeatures;

/*---------------------------------------------------------------------------*\
                           Class topoCellLooper Declaration
\*---------------------------------------------------------------------------*/

class topoCellLooper
:
    public hexCellLooper
{

    // Private Member Functions

        //- In-memory truncate a list
        template<class T>
        static void subsetList
        (
            const label startI,
            const label freeI,
            DynamicList<T>& lst
        );

        //- Walk across superface discarding non-feature points.
        void walkFace
        (
            const cellFeatures& features,
            const label facei,
            const label startEdgeI,
            const label startVertI,
            const label nFeaturePts,

            label& edgeI,
            label& vertI
        ) const;

        //- Returns list of vertices on 'superEdge' i.e. list of edges connected
        // by non-feature points. First and last are feature points, ones
        // in between are not.
        labelList getSuperEdge
        (
            const cellFeatures& features,
            const label facei,
            const label startEdgeI,
            const label startVertI
        ) const;

        // Return non-feature edge from cells' edges that is most
        // perpendicular to refinement direction. Used as starting edge.
        label getAlignedNonFeatureEdge
        (
            const vector& refDir,
            const label celli,
            const cellFeatures& features
        ) const;

        //- Starts from edge and vertex on edge on face (or neighbouring face)
        // and steps either to existing vertex (vertI != -1) or to edge
        // (vertI == -1)
        // by walking point-edge and crossing nFeats featurePoints.
        void walkAcrossFace
        (
            const cellFeatures& features,
            const label facei,
            const label startEdgeI,
            const label startVertI,
            const label nFeats,

            label& edgeI,
            label& vertI
        ) const;

        //- Walks splitcell circumference. Sets loop/loopweights to walk on
        //  outside of cell.
        void walkSplitHex
        (
            const label celli,
            const cellFeatures& features,
            const label fromFacei,
            const label fromEdgeI,
            const label fromVertI,

            DynamicList<label>& loop,
            DynamicList<scalar>& loopWeights
        ) const;


        //- Disallow default bitwise copy construct
        topoCellLooper(const topoCellLooper&);

        //- Disallow default bitwise assignment
        void operator=(const topoCellLooper&);


public:

    //- Runtime type information
    TypeName("topoCellLooper");

    // Static data members

        //- Cos of angle for feature recognition (of splitHexes)
        static const scalar featureCos;


    // Constructors

        //- Construct from components
        topoCellLooper(const polyMesh& mesh);


    //- Destructor
    virtual ~topoCellLooper();


    // Member Functions

        //- Create cut along circumference of celli. Gets current mesh cuts.
        //  Cut along circumference is expressed as loop of cuts plus weights
        //  for cuts along edges (only valid for edge cuts).
        //  Return true if successful cut.
        virtual bool cut
        (
            const vector& refDir,
            const label celli,
            const boolList& vertIsCut,
            const boolList& edgeIsCut,
            const scalarField& edgeWeight,

            labelList& loop,
            scalarField& loopWeights
        ) const;

        //- Same but now also base point of cut provided (instead of always
        //  cell centre)
        virtual bool cut
        (
            const plane& cutPlane,
            const label celli,
            const boolList& vertIsCut,
            const boolList& edgeIsCut,
            const scalarField& edgeWeight,

            labelList& loop,
            scalarField& loopWeights
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
