/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    surfaceRedistributePar

Description
    (Re)distribution of triSurface. Either takes an undecomposed surface
    or an already decomposed surface and redistributes it so that each
    processor has all triangles that overlap its mesh.

Note
    - best decomposition option is hierarchGeomDecomp since
    guarantees square decompositions.
    - triangles might be present on multiple processors.
    - merging uses geometric tolerance so take care with writing precision.

\*---------------------------------------------------------------------------*/

#include "treeBoundBox.H"
#include "FixedList.H"
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "distributedTriSurfaceMesh.H"
#include "mapDistribute.H"
#include "triSurfaceFields.H"
#include "Pair.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Print on master all the per-processor surface stats.
void writeProcStats
(
    const triSurface& s,
    const List<List<treeBoundBox> >& meshBb
)
{
    // Determine surface bounding boxes, faces, points
    List<treeBoundBox> surfBb(Pstream::nProcs());
    {
        surfBb[Pstream::myProcNo()] = treeBoundBox(s.points());
        Pstream::gatherList(surfBb);
        Pstream::scatterList(surfBb);
    }

    labelList nPoints(Pstream::nProcs());
    nPoints[Pstream::myProcNo()] = s.points().size();
    Pstream::gatherList(nPoints);
    Pstream::scatterList(nPoints);

    labelList nFaces(Pstream::nProcs());
    nFaces[Pstream::myProcNo()] = s.size();
    Pstream::gatherList(nFaces);
    Pstream::scatterList(nFaces);

    forAll(surfBb, procI)
    {
        const List<treeBoundBox>& bbs = meshBb[procI];

        Info<< "processor" << procI << nl
            << "\tMesh bounds          : " << bbs[0] << nl;
        for (label i = 1; i < bbs.size(); i++)
        {
            Info<< "\t                       " << bbs[i]<< nl;
        }
        Info<< "\tSurface bounding box : " << surfBb[procI] << nl
            << "\tTriangles            : " << nFaces[procI] << nl
            << "\tVertices             : " << nPoints[procI]
            << endl;
    }
    Info<< endl;
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "redistribute a triSurface"
    );

    argList::validArgs.append("triSurfaceMesh");
    argList::validArgs.append("distributionType");
    argList::addBoolOption
    (
        "keepNonMapped",
        "preserve surface outside of mesh bounds"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();

    const fileName surfFileName = args[1];
    const word distType = args[2];

    Info<< "Reading surface from " << surfFileName << nl << nl
        << "Using distribution method "
        << distributedTriSurfaceMesh::distributionTypeNames_[distType]
        << " " << distType << nl << endl;

    const bool keepNonMapped = args.options().found("keepNonMapped");

    if (keepNonMapped)
    {
        Info<< "Preserving surface outside of mesh bounds." << nl << endl;
    }
    else
    {
        Info<< "Removing surface outside of mesh bounds." << nl << endl;
    }


    if (!Pstream::parRun())
    {
        FatalErrorIn(args.executable())
            << "Please run this program on the decomposed case."
            << " It will read surface " << surfFileName
            << " and decompose it such that it overlaps the mesh bounding box."
            << exit(FatalError);
    }


    #include "createPolyMesh.H"

    Random rndGen(653213);

    // Determine mesh bounding boxes:
    List<List<treeBoundBox> > meshBb(Pstream::nProcs());
    {
        meshBb[Pstream::myProcNo()] = List<treeBoundBox>
        (
            1,
            treeBoundBox
            (
                boundBox(mesh.points(), false)
            ).extend(rndGen, 1e-3)
        );
        Pstream::gatherList(meshBb);
        Pstream::scatterList(meshBb);
    }

    IOobject io
    (
        surfFileName,         // name
        //runTime.findInstance("triSurface", surfFileName),   // instance
        runTime.constant(),   // instance
        "triSurface",         // local
        runTime,              // registry
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    const fileName actualPath(io.filePath());
    fileName localPath(actualPath);
    localPath.replace(runTime.rootPath() + '/', "");

    if (actualPath == io.objectPath())
    {
        Info<< "Loading local (decomposed) surface " << localPath << nl <<endl;
    }
    else
    {
        Info<< "Loading undecomposed surface " << localPath << nl << endl;
    }


    // Create dummy dictionary for bounding boxes if does not exist.
    if (!isFile(actualPath / "Dict"))
    {
        dictionary dict;
        dict.add("bounds", meshBb[Pstream::myProcNo()]);
        dict.add("distributionType", distType);
        dict.add("mergeDistance", SMALL);

        IOdictionary ioDict
        (
            IOobject
            (
                io.name() + "Dict",
                io.instance(),
                io.local(),
                io.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            dict
        );

        Info<< "Writing dummy bounds dictionary to " << ioDict.name()
            << nl << endl;

        ioDict.regIOobject::writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            ioDict.time().writeCompression()
        );
    }


    // Load surface
    distributedTriSurfaceMesh surfMesh(io);
    Info<< "Loaded surface" << nl << endl;


    // Generate a test field
    {
        const triSurface& s = static_cast<const triSurface&>(surfMesh);

        autoPtr<triSurfaceVectorField> fcPtr
        (
            new triSurfaceVectorField
            (
                IOobject
                (
                    surfMesh.searchableSurface::name(),     // name
                    surfMesh.searchableSurface::instance(), // instance
                    surfMesh.searchableSurface::local(),    // local
                    surfMesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                surfMesh,
                dimLength
            )
        );
        triSurfaceVectorField& fc = fcPtr();

        forAll(fc, triI)
        {
            fc[triI] = s[triI].centre(s.points());
        }

        // Steal pointer and store object on surfMesh
        fcPtr.ptr()->store();
    }


    // Write per-processor stats
    Info<< "Before redistribution:" << endl;
    writeProcStats(surfMesh, meshBb);


    // Do redistribution
    Info<< "Redistributing surface" << nl << endl;
    autoPtr<mapDistribute> faceMap;
    autoPtr<mapDistribute> pointMap;
    surfMesh.distribute
    (
        meshBb[Pstream::myProcNo()],
        keepNonMapped,
        faceMap,
        pointMap
    );
    faceMap.clear();
    pointMap.clear();

    Info<< endl;


    // Write per-processor stats
    Info<< "After redistribution:" << endl;
    writeProcStats(surfMesh, meshBb);


    Info<< "Writing surface." << nl << endl;
    surfMesh.searchableSurface::write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
