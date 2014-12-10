/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "searchablePlateFeatures.H"
#include "addToRunTimeSelectionTable.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(searchablePlateFeatures, 0);
addToRunTimeSelectionTable
(
    searchableSurfaceFeatures,
    searchablePlateFeatures,
    dict
);

//! \cond - skip documentation : local scope only
const Foam::label edgesArray[4][2] =
{
    {0, 1}, // 0
    {0, 3},
    {2, 1}, // 2
    {2, 3}
};
//! \endcond

const edgeList searchablePlateFeatures::edges(calcEdges(edgesArray));

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::edgeList Foam::searchablePlateFeatures::calcEdges
(
    const label edgesArray[4][2]
)
{
    edgeList edges(4);
    forAll(edges, edgeI)
    {
        edges[edgeI][0] = edgesArray[edgeI][0];
        edges[edgeI][1] = edgesArray[edgeI][1];
    }
    return edges;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchablePlateFeatures::searchablePlateFeatures
(
    const searchableSurface& surface,
    const dictionary& dict
)
:
    searchableSurfaceFeatures(surface, dict),
    mode_
    (
        extendedFeatureEdgeMesh::sideVolumeTypeNames_
        [
            dict.lookupOrDefault<word>("meshableSide", "inside")
        ]
    )
{
    Info<< indent
        << "    Meshable region = "
        << extendedFeatureEdgeMesh::sideVolumeTypeNames_[mode_]
        << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchablePlateFeatures::~searchablePlateFeatures()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::extendedFeatureEdgeMesh>
Foam::searchablePlateFeatures::features() const
{
    autoPtr<extendedFeatureEdgeMesh> features;

    vectorField faceNormals(1);
    surface().getNormal(List<pointIndexHit>(1), faceNormals);

    vectorField edgeDirections(4);
    labelListList normalDirections(4);

    labelListList edgeNormals(4);
    forAll(edgeNormals, eI)
    {
        edgeNormals[eI].setSize(2, 0);
    }
    edgeNormals[0][0] = 0; edgeNormals[0][1] = 0;
    edgeNormals[1][0] = 0; edgeNormals[1][1] = 0;
    edgeNormals[2][0] = 0; edgeNormals[2][1] = 0;
    edgeNormals[3][0] = 0; edgeNormals[3][1] = 0;

    forAll(edgeDirections, eI)
    {
        edgeDirections[eI] =
            surface().points()()[edges[eI].end()]
          - surface().points()()[edges[eI].start()];

        normalDirections[eI] = labelList(2, 0);
        for (label j = 0; j < 2; ++j)
        {
            const vector cross =
                (faceNormals[edgeNormals[eI][j]] ^ edgeDirections[eI]);
            const vector fC0tofE0 =
                0.5*(max(surface().points()() + min(surface().points()())))
              - surface().points()()[edges[eI].start()];

            normalDirections[eI][j] =
                (
                    (
                        (cross/(mag(cross) + VSMALL))
                      & (fC0tofE0/(mag(fC0tofE0)+ VSMALL))
                    )
                  > 0.0
                    ? 1
                    : -1
                );
        }
    }

    labelListList featurePointNormals(4);
    labelListList featurePointEdges(4);
    forAll(featurePointNormals, pI)
    {
        labelList& ftPtEdges = featurePointEdges[pI];
        ftPtEdges.setSize(2, 0);

        label edgeI = 0;
        forAll(edges, eI)
        {
            const edge& e = edges[eI];

            if (e.start() == pI)
            {
                ftPtEdges[edgeI++] = eI;
            }
            else if (e.end() == pI)
            {
                ftPtEdges[edgeI++] = eI;
            }
        }

        labelList& ftPtNormals = featurePointNormals[pI];
        ftPtNormals.setSize(1, 0);

        ftPtNormals[0] = edgeNormals[ftPtEdges[0]][0];
    }

    labelList regionEdges;

    features.set
    (
        new extendedFeatureEdgeMesh
        (
            IOobject
            (
                surface().name(),
                surface().instance(),
                "extendedFeatureEdgeMesh",
                surface().db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            surface().points(),
            edges,
            4, 4, 4,
            0, 0, 4, 4, // 4 flat edges
            faceNormals,
            List<extendedFeatureEdgeMesh::sideVolumeType>(4, mode_),
            edgeDirections,
            normalDirections,
            edgeNormals,
            featurePointNormals,
            featurePointEdges,
            regionEdges
        )
    );

    return features;
}


// ************************************************************************* //
