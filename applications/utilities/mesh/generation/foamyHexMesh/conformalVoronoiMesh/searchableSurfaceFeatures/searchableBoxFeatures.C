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

#include "searchableBoxFeatures.H"
#include "addToRunTimeSelectionTable.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(searchableBoxFeatures, 0);
addToRunTimeSelectionTable
(
    searchableSurfaceFeatures,
    searchableBoxFeatures,
    dict
);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableBoxFeatures::searchableBoxFeatures
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

Foam::searchableBoxFeatures::~searchableBoxFeatures()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::extendedFeatureEdgeMesh>
Foam::searchableBoxFeatures::features() const
{
    autoPtr<extendedFeatureEdgeMesh> features;

    List<vector> faceNormalsList(treeBoundBox::faceNormals);
    vectorField faceNormals(faceNormalsList);

    vectorField edgeDirections(12);
    labelListList normalDirections(12);

    labelListList edgeNormals(12);
    forAll(edgeNormals, eI)
    {
        edgeNormals[eI].setSize(2, 0);
    }
    edgeNormals[0][0] = 2; edgeNormals[0][1] = 4;
    edgeNormals[1][0] = 1; edgeNormals[1][1] = 4;
    edgeNormals[2][0] = 3; edgeNormals[2][1] = 4;
    edgeNormals[3][0] = 0; edgeNormals[3][1] = 4;
    edgeNormals[4][0] = 2; edgeNormals[4][1] = 5;
    edgeNormals[5][0] = 1; edgeNormals[5][1] = 5;
    edgeNormals[6][0] = 3; edgeNormals[6][1] = 5;
    edgeNormals[7][0] = 0; edgeNormals[7][1] = 5;
    edgeNormals[8][0] = 0; edgeNormals[8][1] = 2;
    edgeNormals[9][0] = 2; edgeNormals[9][1] = 1;
    edgeNormals[10][0] = 1; edgeNormals[10][1] = 3;
    edgeNormals[11][0] = 3; edgeNormals[11][1] = 0;

    tmp<pointField> surfacePointsTmp(surface().points());
    pointField& surfacePoints = surfacePointsTmp();

    forAll(edgeDirections, eI)
    {
        edgeDirections[eI] =
            surfacePoints[treeBoundBox::edges[eI].end()]
          - surfacePoints[treeBoundBox::edges[eI].start()];

        normalDirections[eI] = labelList(2, 0);
        for (label j = 0; j < 2; ++j)
        {
            const vector cross =
                (faceNormals[edgeNormals[eI][j]] ^ edgeDirections[eI]);
            const vector fC0tofE0 =
                0.5*(max(surfacePoints + min(surfacePoints)))
              - surfacePoints[treeBoundBox::edges[eI].start()];

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

    labelListList featurePointNormals(8);
    labelListList featurePointEdges(8);
    forAll(featurePointNormals, pI)
    {
        labelList& ftPtEdges = featurePointEdges[pI];
        ftPtEdges.setSize(3, 0);

        label edgeI = 0;
        forAll(treeBoundBox::edges, eI)
        {
            const edge& e = treeBoundBox::edges[eI];

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
        ftPtNormals.setSize(3, 0);

        ftPtNormals[0] = edgeNormals[ftPtEdges[0]][0];
        ftPtNormals[1] = edgeNormals[ftPtEdges[0]][1];
        ftPtNormals[2] = edgeNormals[ftPtEdges[1]][0];
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
            treeBoundBox::edges,
            8, 8, 8,
            12, 12, 12, 12,
            faceNormals,
            List<extendedFeatureEdgeMesh::sideVolumeType>(12, mode_),
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
