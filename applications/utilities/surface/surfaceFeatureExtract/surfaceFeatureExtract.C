/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

\*---------------------------------------------------------------------------*/

#include "surfaceFeatureExtract.H"
#include "argList.H"
#include "Time.H"
#include "featureEdgeMesh.H"
#include "vtkSurfaceWriter.H"
#include "IOdictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "extract and write surface features to file"
    );
    argList::noParallel();

    #include "addDictOption.H"

    #include "setRootCase.H"
    #include "createTime.H"

    const word dictName("surfaceFeatureExtractDict");
    #include "setSystemRunTimeDictionaryIO.H"

    Info<< "Reading " << dictName << nl << endl;

    const IOdictionary dict(dictIO);

    forAllConstIter(dictionary, dict, iter)
    {
        if (!iter().isDict())
        {
            continue;
        }

        const dictionary& surfaceDict = iter().dict();

        if (!surfaceDict.found("extractionMethod"))
        {
            continue;
        }

        const word extractionMethod = surfaceDict.lookup("extractionMethod");

        const fileName surfFileName = iter().keyword();
        const fileName sFeatFileName = surfFileName.lessExt().name();

        Info<< "Surface            : " << surfFileName << nl << endl;

        const Switch writeVTK =
            surfaceDict.lookupOrDefault<Switch>("writeVTK", "off");
        const Switch writeObj =
            surfaceDict.lookupOrDefault<Switch>("writeObj", "off");

        const Switch curvature =
            surfaceDict.lookupOrDefault<Switch>("curvature", "off");
        const Switch featureProximity =
            surfaceDict.lookupOrDefault<Switch>("featureProximity", "off");
        const Switch closeness =
            surfaceDict.lookupOrDefault<Switch>("closeness", "off");


        Info<< nl << "Feature line extraction is only valid on closed manifold "
            << "surfaces." << endl;

        // Read
        // ~~~~

        triSurface surf(runTime.constantPath()/"triSurface"/surfFileName);

        Info<< "Statistics:" << endl;
        surf.writeStats(Info);
        Info<< endl;

        const faceList faces(surf.faces());

        // Either construct features from surface & featureAngle or read set.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        autoPtr<surfaceFeatures> set;

        scalar includedAngle = 0.0;

        if (extractionMethod == "extractFromFile")
        {
            const dictionary& extractFromFileDict =
                surfaceDict.subDict("extractFromFileCoeffs");

            const fileName featureEdgeFile =
                extractFromFileDict.lookup("featureEdgeFile");

            const Switch geometricTestOnly =
                extractFromFileDict.lookupOrDefault<Switch>
                (
                    "geometricTestOnly",
                    "no"
                );

            edgeMesh eMesh(featureEdgeFile);

            // Sometimes duplicate edges are present. Remove them.
            eMesh.mergeEdges();

            Info<< nl << "Reading existing feature edges from file "
                << featureEdgeFile << nl
                << "Selecting edges purely based on geometric tests: "
                << geometricTestOnly.asText() << endl;

            set.set
            (
                new surfaceFeatures
                (
                    surf,
                    eMesh.points(),
                    eMesh.edges(),
                    1e-6,
                    geometricTestOnly
                )
            );
        }
        else if (extractionMethod == "extractFromSurface")
        {
            const dictionary& extractFromSurfaceDict =
                surfaceDict.subDict("extractFromSurfaceCoeffs");

            includedAngle =
                readScalar(extractFromSurfaceDict.lookup("includedAngle"));

            const Switch geometricTestOnly =
                extractFromSurfaceDict.lookupOrDefault<Switch>
                (
                    "geometricTestOnly",
                    "no"
                );

            Info<< nl
                << "Constructing feature set from included angle "
                << includedAngle << nl
                << "Selecting edges purely based on geometric tests: "
                << geometricTestOnly.asText() << endl;

            set.set
            (
                new surfaceFeatures
                (
                    surf,
                    includedAngle,
                    0,
                    0,
                    geometricTestOnly
                )
            );
        }
        else
        {
            FatalErrorInFunction
                << "No initial feature set. Provide either one"
                << " of extractFromFile (to read existing set)" << nl
                << " or extractFromSurface (to construct new set from angle)"
                << exit(FatalError);
        }


        // Trim set
        // ~~~~~~~~

        if (surfaceDict.isDict("trimFeatures"))
        {
            dictionary trimDict = surfaceDict.subDict("trimFeatures");

            scalar minLen =
                trimDict.lookupOrAddDefault<scalar>("minLen", -great);

            label minElem = trimDict.lookupOrAddDefault<label>("minElem", 0);

            // Trim away small groups of features
            if (minElem > 0 || minLen > 0)
            {
                Info<< "Removing features of length < "
                    << minLen << endl;
                Info<< "Removing features with number of edges < "
                    << minElem << endl;

                set().trimFeatures(minLen, minElem, includedAngle);
            }
        }


        // Subset
        // ~~~~~~

        // Convert to marked edges, points
        List<surfaceFeatures::edgeStatus> edgeStat(set().toStatus());

        if (surfaceDict.isDict("subsetFeatures"))
        {
            const dictionary& subsetDict = surfaceDict.subDict
            (
                "subsetFeatures"
            );

            if (subsetDict.found("insideBox"))
            {
                treeBoundBox bb(subsetDict.lookup("insideBox")());

                Info<< "Removing all edges outside bb " << bb << endl;
                dumpBox(bb, "subsetBox.obj");

                deleteBox(surf, bb, false, edgeStat);
            }
            else if (subsetDict.found("outsideBox"))
            {
                treeBoundBox bb(subsetDict.lookup("outsideBox")());

                Info<< "Removing all edges inside bb " << bb << endl;
                dumpBox(bb, "deleteBox.obj");

                deleteBox(surf, bb, true, edgeStat);
            }

            const Switch nonManifoldEdges =
                subsetDict.lookupOrDefault<Switch>("nonManifoldEdges", "yes");

            if (!nonManifoldEdges)
            {
                Info<< "Removing all non-manifold edges"
                    << " (edges with > 2 connected faces) unless they"
                    << " cross multiple regions" << endl;

                forAll(edgeStat, edgeI)
                {
                    const labelList& eFaces = surf.edgeFaces()[edgeI];

                    if
                    (
                        eFaces.size() > 2
                     && edgeStat[edgeI] == surfaceFeatures::REGION
                     && (eFaces.size() % 2) == 0
                    )
                    {
                        edgeStat[edgeI] = checkFlatRegionEdge
                        (
                            surf,
                            1e-5,   //tol,
                            includedAngle,
                            edgeI
                        );
                    }
                }
            }

            const Switch openEdges =
                subsetDict.lookupOrDefault<Switch>("openEdges", "yes");

            if (!openEdges)
            {
                Info<< "Removing all open edges"
                    << " (edges with 1 connected face)" << endl;

                forAll(edgeStat, edgeI)
                {
                    if (surf.edgeFaces()[edgeI].size() == 1)
                    {
                        edgeStat[edgeI] = surfaceFeatures::NONE;
                    }
                }
            }

            if (subsetDict.found("plane"))
            {
                plane cutPlane(subsetDict.lookup("plane")());

                deleteEdges(surf, cutPlane, edgeStat);

                Info<< "Only edges that intersect the plane with normal "
                    << cutPlane.normal()
                    << " and base point " << cutPlane.refPoint()
                    << " will be included as feature edges."<< endl;
            }
        }


        surfaceFeatures newSet(surf);
        newSet.setFromStatus(edgeStat, includedAngle);

        Info<< nl
            << "Initial feature set:" << nl
            << "    feature points : " << newSet.featurePoints().size() << nl
            << "    feature edges  : " << newSet.featureEdges().size() << nl
            << "    of which" << nl
            << "        region edges   : " << newSet.nRegionEdges() << nl
            << "        external edges : " << newSet.nExternalEdges() << nl
            << "        internal edges : " << newSet.nInternalEdges() << nl
            << endl;

        boolList surfBaffleRegions(surf.patches().size(), false);

        wordList surfBaffleNames;
        surfaceDict.readIfPresent("baffles", surfBaffleNames);

        forAll(surf.patches(), pI)
        {
            const word& name = surf.patches()[pI].name();

            if (findIndex(surfBaffleNames, name) != -1)
            {
                Info<< "Adding baffle region " << name << endl;
                surfBaffleRegions[pI] = true;
            }
        }

        // Extracting and writing a extendedFeatureEdgeMesh
        extendedFeatureEdgeMesh feMesh
        (
            newSet,
            runTime,
            sFeatFileName + ".extendedFeatureEdgeMesh",
            surfBaffleRegions
        );


        if (surfaceDict.isDict("addFeatures"))
        {
            const word addFeName = surfaceDict.subDict("addFeatures")["name"];
            Info<< "Adding (without merging) features from " << addFeName
                << nl << endl;

            extendedFeatureEdgeMesh addFeMesh
            (
                IOobject
                (
                    addFeName,
                    runTime.time().constant(),
                    "extendedFeatureEdgeMesh",
                    runTime.time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            Info<< "Read " << addFeMesh.name() << nl;
            writeStats(addFeMesh, Info);

            feMesh.add(addFeMesh);
        }


        Info<< nl
            << "Final feature set:" << nl;
        writeStats(feMesh, Info);

        Info<< nl << "Writing extendedFeatureEdgeMesh to "
            << feMesh.objectPath() << endl;

        mkDir(feMesh.path());

        if (writeObj)
        {
            feMesh.writeObj(feMesh.path()/surfFileName.lessExt().name());
        }

        feMesh.write();

        // Write a featureEdgeMesh for backwards compatibility
        featureEdgeMesh bfeMesh
        (
            IOobject
            (
                surfFileName.lessExt().name() + ".eMesh",   // name
                runTime.constant(),                         // instance
                "triSurface",
                runTime,                                    // registry
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            feMesh.points(),
            feMesh.edges()
        );

        Info<< nl << "Writing featureEdgeMesh to "
            << bfeMesh.objectPath() << endl;

        bfeMesh.regIOobject::write();


        // Find distance between close features
        if (closeness)
        {
            Info<< nl << "Extracting internal and external closeness of "
                << "surface." << endl;

            Info<< "externalToleranceCosAngle: " << externalToleranceCosAngle
                << nl
                << "internalToleranceCosAngle: " << internalToleranceCosAngle
                << endl;

            extractCloseness(sFeatFileName, runTime, surf, writeVTK);
            extractPointCloseness(sFeatFileName, runTime, surf, writeVTK);
        }


        if (curvature)
        {
            Info<< nl << "Extracting curvature of surface at the points."
                << endl;

            const vectorField pointNormals(surf.pointNormals2());
            triadField pointCoordSys = calcVertexCoordSys(surf, pointNormals);

            triSurfacePointScalarField k = calcCurvature
            (
                sFeatFileName,
                runTime,
                surf,
                pointNormals,
                pointCoordSys
            );

            k.write();

            if (writeVTK)
            {
                vtkSurfaceWriter().write
                (
                    runTime.constantPath()/"triSurface",// outputDir
                    sFeatFileName,                      // surfaceName
                    surf.points(),
                    faces,
                    "curvature",                        // fieldName
                    k,
                    true,                               // isNodeValues
                    true                                // verbose
                );
            }
        }


        if (featureProximity)
        {
            Info<< nl << "Extracting proximity of close feature points and "
                << "edges to the surface" << endl;

            const scalar searchDistance =
                readScalar(surfaceDict.lookup("maxFeatureProximity"));

            scalarField featureProximity(surf.size(), searchDistance);

            forAll(surf, fI)
            {
                const triPointRef& tri = surf[fI].tri(surf.points());
                const point& triCentre = tri.circumCentre();

                const scalar radiusSqr = min
                (
                    sqr(4*tri.circumRadius()),
                    sqr(searchDistance)
                );

                List<pointIndexHit> hitList;

                feMesh.allNearestFeatureEdges(triCentre, radiusSqr, hitList);

                featureProximity[fI] =
                    calcProximityOfFeatureEdges
                    (
                        feMesh,
                        hitList,
                        featureProximity[fI]
                    );

                feMesh.allNearestFeaturePoints(triCentre, radiusSqr, hitList);

                featureProximity[fI] =
                    calcProximityOfFeaturePoints
                    (
                        hitList,
                        featureProximity[fI]
                    );
            }

            triSurfaceScalarField featureProximityField
            (
                IOobject
                (
                    sFeatFileName + ".featureProximity",
                    runTime.constant(),
                    "triSurface",
                    runTime,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                surf,
                dimLength,
                featureProximity
            );

            featureProximityField.write();

            if (writeVTK)
            {
                vtkSurfaceWriter().write
                (
                    runTime.constantPath()/"triSurface",// outputDir
                    sFeatFileName,                      // surfaceName
                    surf.points(),
                    faces,
                    "featureProximity",                 // fieldName
                    featureProximity,
                    false,                              // isNodeValues
                    true                                // verbose
                );
            }
        }

        Info<< endl;
    }

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
