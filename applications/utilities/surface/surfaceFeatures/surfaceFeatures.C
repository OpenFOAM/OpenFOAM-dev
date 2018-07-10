/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

Description
    Identifies features in a surface geometry and writes them to file,
    based on control parameters specified by the user.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "triSurfaceMesh.H"
#include "featureEdgeMesh.H"
#include "extendedFeatureEdgeMesh.H"
#include "surfaceFeatures.H"
#include "triSurfaceFields.H"
#include "vtkSurfaceWriter.H"
#include "IOdictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    autoPtr<surfaceFeatures> extractFromFile
    (
        const fileName& featureEdgeFile,
        const triSurface& surf,
        const Switch& geometricTestOnly
    )
    {
        edgeMesh eMesh(featureEdgeFile);

        // Sometimes duplicate edges are present. Remove them.
        eMesh.mergeEdges();

        Info<< nl << "Reading existing feature edges from file "
            << featureEdgeFile << nl
            << "Selecting edges purely based on geometric tests: "
            << geometricTestOnly.asText() << endl;

        return  autoPtr<surfaceFeatures>
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


    autoPtr<surfaceFeatures> extractFromSurface
    (
        const triSurface& surf,
        const Switch& geometricTestOnly,
        const scalar includedAngle
    )
    {
        Info<< nl
            << "Constructing feature set from included angle "
            << includedAngle << nl
            << "Selecting edges purely based on geometric tests: "
            << geometricTestOnly.asText() << endl;

        return  autoPtr<surfaceFeatures>
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


    autoPtr<surfaceFeatures> surfaceFeatureSet
    (
        const fileName& surfaceFileName,
        const triSurface& surf,
        const dictionary& dict,
        const scalar includedAngle
    )
    {
        const Switch geometricTestOnly = dict.lookupOrDefault<Switch>
        (
            "geometricTestOnly",
            "no"
        );

        if (dict.found("files"))
        {
            HashTable<fileName, fileName> fileNames(dict.lookup("files"));

            if (fileNames.found(surfaceFileName))
            {
                return extractFromFile
                (
                    fileNames[surfaceFileName],
                    surf,
                    geometricTestOnly
                );
            }
            else
            {
                return extractFromSurface
                (
                    surf,
                    geometricTestOnly,
                    includedAngle
                );
            }
        }
        else
        {
            return extractFromSurface
            (
                surf,
                geometricTestOnly,
                includedAngle
            );
        }
    }


    void extractFeatures
    (
        const fileName& surfaceFileName,
        const Time& runTime,
        const dictionary& dict
    )
    {
        const fileName sFeatFileName = surfaceFileName.lessExt().name();

        Info<< "Surface            : " << surfaceFileName << nl << endl;

        const Switch writeVTK =
            dict.lookupOrDefault<Switch>("writeVTK", "off");
        const Switch writeObj =
            dict.lookupOrDefault<Switch>("writeObj", "off");
        const Switch verboseObj =
            dict.lookupOrDefault<Switch>("verboseObj", "off");

        const Switch curvature =
            dict.lookupOrDefault<Switch>("curvature", "off");
        const Switch featureProximity =
            dict.lookupOrDefault<Switch>("featureProximity", "off");
        const Switch closeness =
            dict.lookupOrDefault<Switch>("closeness", "off");


        Info<< nl << "Feature line extraction is only valid on closed manifold "
            << "surfaces." << endl;

        // Read
        // ~~~~

        triSurface surf(runTime.constantPath()/"triSurface"/surfaceFileName);

        Info<< "Statistics:" << endl;
        surf.writeStats(Info);
        Info<< endl;

        const faceList faces(surf.faces());

        // Either construct features from surface & featureAngle or read set.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        const scalar includedAngle = readScalar(dict.lookup("includedAngle"));

        autoPtr<surfaceFeatures> set
        (
            surfaceFeatureSet
            (
                surfaceFileName,
                surf,
                dict,
                includedAngle
            )
        );

        // Trim set
        // ~~~~~~~~

        if (dict.isDict("trimFeatures"))
        {
            dictionary trimDict = dict.subDict("trimFeatures");

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

        if (dict.isDict("subsetFeatures"))
        {
            const dictionary& subsetDict = dict.subDict
            (
                "subsetFeatures"
            );

            if (subsetDict.found("insideBox"))
            {
                treeBoundBox bb(subsetDict.lookup("insideBox")());

                Info<< "Selecting edges inside bb " << bb;
                if (writeObj)
                {
                    Info << " see insideBox.obj";
                    bb.writeOBJ("insideBox.obj");
                }
                Info<< endl;

                selectBox(surf, bb, true, edgeStat);
            }
            else if (subsetDict.found("outsideBox"))
            {
                treeBoundBox bb(subsetDict.lookup("outsideBox")());

                Info<< "Removing all edges inside bb " << bb;
                if (writeObj)
                {
                    Info<< " see outsideBox.obj" << endl;
                    bb.writeOBJ("outsideBox.obj");
                }
                Info<< endl;

                selectBox(surf, bb, false, edgeStat);
            }

            const Switch nonManifoldEdges =
                subsetDict.lookupOrDefault<Switch>("nonManifoldEdges", "yes");

            if (!nonManifoldEdges)
            {
                Info<< "Removing all non-manifold edges"
                    << " (edges with > 2 connected faces) unless they"
                    << " cross multiple regions" << endl;

                selectManifoldEdges(surf, 1e-5, includedAngle, edgeStat);
            }

            const Switch openEdges =
                subsetDict.lookupOrDefault<Switch>("openEdges", "yes");

            if (!openEdges)
            {
                Info<< "Removing all open edges"
                    << " (edges with 1 connected face)" << endl;

                forAll(edgeStat, edgei)
                {
                    if (surf.edgeFaces()[edgei].size() == 1)
                    {
                        edgeStat[edgei] = surfaceFeatures::NONE;
                    }
                }
            }

            if (subsetDict.found("plane"))
            {
                const plane cutPlane(subsetDict.lookup("plane")());

                selectCutEdges(surf, cutPlane, edgeStat);

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
        dict.readIfPresent("baffles", surfBaffleNames);

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


        if (dict.isDict("addFeatures"))
        {
            const word addFeName = dict.subDict("addFeatures")["name"];
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
            addFeMesh.writeStats(Info);

            feMesh.add(addFeMesh);
        }


        Info<< nl
            << "Final feature set:" << nl;
        feMesh.writeStats(Info);

        Info<< nl << "Writing extendedFeatureEdgeMesh to "
            << feMesh.objectPath() << endl;

        mkDir(feMesh.path());

        if (writeObj)
        {
            feMesh.writeObj
            (
                feMesh.path()/surfaceFileName.lessExt().name(),
                verboseObj
            );
        }

        feMesh.write();

        // Write a featureEdgeMesh for backwards compatibility
        featureEdgeMesh bfeMesh
        (
            IOobject
            (
                surfaceFileName.lessExt().name() + ".eMesh",
                runTime.constant(),
                "triSurface",
                runTime,
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

            // Searchable triSurface
            const triSurfaceMesh searchSurf
            (
                IOobject
                (
                    sFeatFileName + ".closeness",
                    runTime.constant(),
                    "triSurface",
                    runTime
                ),
                surf
            );

            {
                Pair<tmp<triSurfaceScalarField>> closenessFields
                (
                    searchSurf.extractCloseness()
                );

                closenessFields.first()->write();
                closenessFields.second()->write();

                if (writeVTK)
                {
                    const faceList faces(searchSurf.faces());

                    vtkSurfaceWriter().write
                    (
                        runTime.constantPath()/"triSurface",// outputDir
                        searchSurf.objectRegistry::name(),  // surfaceName
                        searchSurf.points(),
                        faces,
                        "internalCloseness",                // fieldName
                        closenessFields.first(),
                        false,                              // isNodeValues
                        true                                // verbose
                    );

                    vtkSurfaceWriter().write
                    (
                        runTime.constantPath()/"triSurface",// outputDir
                        searchSurf.objectRegistry::name(),  // surfaceName
                        searchSurf.points(),
                        faces,
                        "externalCloseness",                // fieldName
                        closenessFields.second(),
                        false,                              // isNodeValues
                        true                                // verbose
                    );
                }
            }

            {
                Pair<tmp<triSurfacePointScalarField >> closenessFields
                (
                    searchSurf.extractPointCloseness()
                );

                closenessFields.first()->write();
                closenessFields.second()->write();

                if (writeVTK)
                {
                    const faceList faces(searchSurf.faces());
                    const Map<label>& meshPointMap = searchSurf.meshPointMap();

                    const triSurfacePointScalarField&
                        internalClosenessPointField = closenessFields.first();

                    const triSurfacePointScalarField&
                        externalClosenessPointField = closenessFields.second();

                    scalarField internalCloseness(searchSurf.nPoints(), great);
                    scalarField externalCloseness(searchSurf.nPoints(), great);

                    forAll(meshPointMap, pi)
                    {
                        internalCloseness[pi] =
                            internalClosenessPointField[meshPointMap[pi]];

                        externalCloseness[pi] =
                            externalClosenessPointField[meshPointMap[pi]];
                    }

                    vtkSurfaceWriter().write
                    (
                        runTime.constantPath()/"triSurface",// outputDir
                        searchSurf.objectRegistry::name(),  // surfaceName
                        searchSurf.points(),
                        faces,
                        "internalPointCloseness",           // fieldName
                        internalCloseness,
                        true,                               // isNodeValues
                        true                                // verbose
                    );

                    vtkSurfaceWriter().write
                    (
                        runTime.constantPath()/"triSurface",// outputDir
                        searchSurf.objectRegistry::name(),  // surfaceName
                        searchSurf.points(),
                        faces,
                        "externalPointCloseness",           // fieldName
                        externalCloseness,
                        true,                               // isNodeValues
                        true                                // verbose
                    );
                }
            }
        }


        if (curvature)
        {
            Info<< nl << "Extracting curvature of surface at the points."
                << endl;

            triSurfacePointScalarField k
            (
                IOobject
                (
                    sFeatFileName + ".curvature",
                    runTime.constant(),
                    "triSurface",
                    runTime
                ),
                surf,
                dimLength,
                surf.curvature()
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
                readScalar(dict.lookup("maxFeatureProximity"));

            scalarField featureProximity(surf.size(), searchDistance);

            forAll(surf, fi)
            {
                const triPointRef& tri = surf[fi].tri(surf.points());
                const point& triCentre = tri.circumCentre();

                const scalar radiusSqr = min
                (
                    sqr(4*tri.circumRadius()),
                    sqr(searchDistance)
                );

                pointIndexHitList hitList;

                feMesh.allNearestFeatureEdges(triCentre, radiusSqr, hitList);
                featureProximity[fi] = min
                (
                    feMesh.minDisconnectedDist(hitList),
                    featureProximity[fi]
                );

                feMesh.allNearestFeaturePoints(triCentre, radiusSqr, hitList);
                featureProximity[fi] = min
                (
                    minDist(hitList),
                    featureProximity[fi]
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


    void extractFeatures
    (
        const fileNameList& surfaceFileNames,
        const Time& runTime,
        const dictionary& dict
    )
    {
        forAll(surfaceFileNames, i)
        {
            extractFeatures(surfaceFileNames[i], runTime, dict);
        }
    }
}


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

    const word dictName("surfaceFeaturesDict");
    #include "setSystemRunTimeDictionaryIO.H"

    Info<< "Reading " << dictName << nl << endl;

    const IOdictionary dict(dictIO);

    if (dict.found("surfaces"))
    {
        extractFeatures
        (
            fileNameList(dict.lookup("surfaces")),
            runTime,
            dict
        );
    }
    else
    {
        forAllConstIter(dictionary, dict, iter)
        {
            if (!iter().isDict())
            {
                continue;
            }

            extractFeatures
            (
                fileNameList(iter().dict().lookup("surfaces")),
                runTime,
                iter().dict()
            );
        }
    }

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}



// ************************************************************************* //
