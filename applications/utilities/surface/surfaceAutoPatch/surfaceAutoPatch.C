/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

Application
    surfaceAutoPatch

Description
    Patches surface according to feature angle. Like autoPatch.

\*---------------------------------------------------------------------------*/

#include "triangle.H"
#include "triSurface.H"
#include "argList.H"
#include "surfaceFeatures.H"
#include "treeBoundBox.H"
#include "meshTools.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    argList::validArgs.append("surface file");
    argList::validArgs.append("output surface file");
    argList::validArgs.append("includedAngle [0..180]");
    argList args(argc, argv);

    const fileName inFileName  = args[1];
    const fileName outFileName = args[2];
    const scalar includedAngle = degToRad(args.argRead<scalar>(3));

    Info<< "Surface        : " << inFileName << nl << endl;


    // Read
    // ~~~~

    Info<< "Reading : " << inFileName << endl;
    triSurface surf(inFileName);

    Info<< "Read surface:" << endl;
    surf.writeStats(Info);
    Info<< endl;



    // Construct features from surface&featureangle
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Constructing feature set from included angle " << includedAngle
        << endl;

    surfaceFeatures set(surf, includedAngle);

    Info<< nl
        << "Feature set:" << nl
        << "    feature points : " << set.featurePoints().size() << nl
        << "    feature edges  : " << set.featureEdges().size() << nl
        << "    of which" << nl
        << "        region edges   : " << set.nRegionEdges() << nl
        << "        external edges : " << set.nExternalEdges() << nl
        << "        internal edges : " << set.nInternalEdges() << nl
        << endl;

    // Get per-edge status.
    boolList borderEdge(surf.nEdges(), false);
    forAll(set.featureEdges(), i)
    {
        borderEdge[set.featureEdges()[i]] = true;
    }

    labelList faceRegion(surf.size());
    label nRegions = surf.markZones(borderEdge, faceRegion);

    // Reregion triangles.
    forAll(surf, i)
    {
        surf[i].region() = faceRegion[i];
    }

    // Create some patches
    surf.patches().setSize(nRegions);

    forAll(surf.patches(), patchi)
    {
        surf.patches()[patchi].name() = "patch" + Foam::name(patchi);
        surf.patches()[patchi].geometricType() = "empty";
    }


    Info<< "Writing : " << outFileName << endl;
    surf.write(outFileName, true);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
