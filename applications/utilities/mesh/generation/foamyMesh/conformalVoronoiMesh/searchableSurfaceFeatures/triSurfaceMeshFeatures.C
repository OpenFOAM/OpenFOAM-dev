/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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

#include "triSurfaceMeshFeatures.H"
#include "addToRunTimeSelectionTable.H"
#include "triSurfaceMesh.H"
#include "surfaceFeatures.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(triSurfaceMeshFeatures, 0);
addToRunTimeSelectionTable
(
    searchableSurfaceFeatures,
    triSurfaceMeshFeatures,
    dict
);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceMeshFeatures::triSurfaceMeshFeatures
(
    const searchableSurface& surface,
    const dictionary& dict
)
:
    searchableSurfaceFeatures(surface, dict),
    includedAngle_(dict.lookup<scalar>("includedAngle")),
    mode_
    (
        extendedFeatureEdgeMesh::sideVolumeTypeNames_
        [
            dict.lookupOrDefault<word>("meshableSide", "inside")
        ]
    )
{
    Info<< indent
        << "    Included angle = " << includedAngle_ << nl
        << "    Meshable region = "
        << extendedFeatureEdgeMesh::sideVolumeTypeNames_[mode_]
        << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceMeshFeatures::~triSurfaceMeshFeatures()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::extendedFeatureEdgeMesh>
Foam::triSurfaceMeshFeatures::features() const
{
    autoPtr<extendedFeatureEdgeMesh> features;

    const triSurfaceMesh& surfMesh = refCast<const triSurfaceMesh>(surface());

    surfaceFeatures sFeat(surfMesh, includedAngle_);

    // TODO: Need to read on a per region basis
    boolList surfBaffleRegions
    (
        surfMesh.patches().size(),
        (mode_ == extendedFeatureEdgeMesh::BOTH ? true : false)
    );

    features.set
    (
        new extendedFeatureEdgeMesh
        (
            sFeat,
            surface().db(),
            surface().name() + ".extendedFeatureEdgeMesh",
            surfBaffleRegions
        )
    );

    return features;
}


// ************************************************************************* //
