/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "mergePoints.H"
#include "Time.H"
#include "ListOps.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "surfMesh.H"
#include "primitivePatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Face>
inline bool Foam::MeshedSurface<Face>::isTri()
{
    return false;
}


template<class Face>
Foam::wordHashSet Foam::MeshedSurface<Face>::readTypes()
{
    return wordHashSet(*fileExtensionConstructorTablePtr_);
}


template<class Face>
Foam::wordHashSet Foam::MeshedSurface<Face>::writeTypes()
{
    return wordHashSet(*writefileExtensionMemberFunctionTablePtr_);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Face>
bool Foam::MeshedSurface<Face>::canReadType
(
    const word& ext,
    const bool verbose
)
{
    return fileFormats::surfaceFormatsCore::checkSupport
    (
        readTypes() | FriendType::readTypes(),
        ext,
        verbose,
        "reading"
   );
}


template<class Face>
bool Foam::MeshedSurface<Face>::canWriteType
(
    const word& ext,
    const bool verbose
)
{
    return fileFormats::surfaceFormatsCore::checkSupport
    (
        writeTypes() | ProxyType::writeTypes(),
        ext,
        verbose,
        "writing"
    );
}


template<class Face>
bool Foam::MeshedSurface<Face>::canRead
(
    const fileName& name,
    const bool verbose
)
{
    word ext = name.ext();
    if (ext == "gz")
    {
        ext = name.lessExt().ext();
    }
    return canReadType(ext, verbose);
}


template<class Face>
void Foam::MeshedSurface<Face>::write
(
    const fileName& name,
    const MeshedSurface<Face>& surf
)
{
    if (debug)
    {
        Info<< "MeshedSurface::write"
            "(const fileName&, const MeshedSurface&) : "
            "writing to " << name
            << endl;
    }

    const word ext = name.ext();

    typename writefileExtensionMemberFunctionTable::iterator mfIter =
        writefileExtensionMemberFunctionTablePtr_->find(ext);

    if (mfIter == writefileExtensionMemberFunctionTablePtr_->end())
    {
        // no direct writer, delegate to proxy if possible
        wordHashSet supported = ProxyType::writeTypes();

        if (supported.found(ext))
        {
            MeshedSurfaceProxy<Face>(surf).write(name);
        }
        else
        {
            FatalErrorIn
            (
                "MeshedSurface::write"
                "(const fileName&, const MeshedSurface&)"
            )   << "Unknown file extension " << ext << nl << nl
                << "Valid types are :" << endl
                << (supported | writeTypes())
                << exit(FatalError);
        }
    }
    else
    {
        mfIter()(name, surf);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface()
:
    ParentType(List<Face>(), pointField())
{}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face> >& faceLst,
    const Xfer<surfZoneList>& zoneLst
)
:
    ParentType(List<Face>(), pointField()),
    zones_()
{
    reset(pointLst, faceLst, zoneLst);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face> >& faceLst,
    const labelUList& zoneSizes,
    const UList<word>& zoneNames
)
:
    ParentType(List<Face>(), pointField())
{
    reset(pointLst, faceLst, Xfer<surfZoneList>());

    if (zoneSizes.size())
    {
        if (zoneNames.size())
        {
            addZones(zoneSizes, zoneNames);
        }
        else
        {
            addZones(zoneSizes);
        }
    }
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const MeshedSurface<Face>& surf
)
:
    ParentType(surf.faces(), surf.points()),
    zones_(surf.surfZones())
{}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const UnsortedMeshedSurface<Face>& surf
)
:
    ParentType(List<Face>(), surf.points())
{
    labelList faceMap;
    this->storedZones() = surf.sortedZones(faceMap);

    const List<Face>& origFaces = surf.faces();
    List<Face> newFaces(origFaces.size());

    forAll(newFaces, faceI)
    {
        newFaces[faceMap[faceI]] = origFaces[faceI];
    }

    this->storedFaces().transfer(newFaces);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(const surfMesh& mesh)
:
    ParentType(List<Face>(), pointField())
{
    // same face type as surfMesh
    MeshedSurface<face> surf
    (
        xferCopy(mesh.points()),
        xferCopy(mesh.faces()),
        xferCopy(mesh.surfZones())
    );

    this->transcribe(surf);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const polyBoundaryMesh& bMesh,
    const bool useGlobalPoints
)
:
    ParentType(List<Face>(), pointField())
{
    const polyMesh& mesh = bMesh.mesh();
    const polyPatchList& bPatches = bMesh;

    // Get a single patch for all boundaries
    primitivePatch allBoundary
    (
        SubList<face>
        (
            mesh.faces(),
            mesh.nFaces() - mesh.nInternalFaces(),
            mesh.nInternalFaces()
        ),
        mesh.points()
    );

    // use global/local points:
    const pointField& bPoints =
    (
        useGlobalPoints ? mesh.points() : allBoundary.localPoints()
    );

    // global/local face addressing:
    const List<Face>& bFaces =
    (
        useGlobalPoints ? allBoundary : allBoundary.localFaces()
    );


    // create zone list
    surfZoneList newZones(bPatches.size());

    label startFaceI = 0;
    label nZone = 0;
    forAll(bPatches, patchI)
    {
        const polyPatch& p = bPatches[patchI];

        if (p.size())
        {
            newZones[nZone] = surfZone
            (
                p.name(),
                p.size(),
                startFaceI,
                nZone
            );

            nZone++;
            startFaceI += p.size();
        }
    }

    newZones.setSize(nZone);

    // same face type as the polyBoundaryMesh
    MeshedSurface<face> surf
    (
        xferCopy(bPoints),
        xferCopy(bFaces),
        xferMove(newZones)
    );

    this->transcribe(surf);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const fileName& name,
    const word& ext
)
:
    ParentType(List<Face>(), pointField())
{
    read(name, ext);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(const fileName& name)
:
    ParentType(List<Face>(), pointField())
{
    read(name);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const Time& t,
    const word& surfName
)
:
    ParentType(List<Face>(), pointField())
{
    surfMesh mesh
    (
        IOobject
        (
            "dummyName",
            t.timeName(),
            t,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        ),
        surfName
    );

    // same face type as surfMesh
    MeshedSurface<face> surf
    (
        xferMove(mesh.storedPoints()),
        xferMove(mesh.storedFaces()),
        xferMove(mesh.storedZones())
    );

    this->transcribe(surf);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const Xfer<UnsortedMeshedSurface<Face> >& surf
)
:
    ParentType(List<Face>(), pointField())
{
    transfer(surf());
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const Xfer<MeshedSurface<Face> >& surf
)
:
    ParentType(List<Face>(), pointField())
{
    transfer(surf());
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurface<Face>::~MeshedSurface()
{}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::remapFaces
(
    const labelUList& faceMap
)
{
    // recalculate the zone start/size
    if (&faceMap && faceMap.size())
    {
        surfZoneList& zones = storedZones();

        if (zones.size() == 1)
        {
            // optimized for single zone case
            zones[0].size() = faceMap.size();
        }
        else if (zones.size())
        {
            label newFaceI = 0;
            label origEndI = 0;
            forAll(zones, zoneI)
            {
                surfZone& zone = zones[zoneI];

                // adjust zone start
                zone.start() = newFaceI;
                origEndI += zone.size();

                for (label faceI = newFaceI; faceI < faceMap.size(); ++faceI)
                {
                    if (faceMap[faceI] < origEndI)
                    {
                        ++newFaceI;
                    }
                    else
                    {
                        break;
                    }
                }

                // adjust zone size
                zone.size() = newFaceI - zone.start();
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::clear()
{
    ParentType::clearOut();

    storedPoints().clear();
    storedFaces().clear();
    storedZones().clear();
}


template<class Face>
void Foam::MeshedSurface<Face>::movePoints(const pointField& newPoints)
{
    // Adapt for new point position
    ParentType::movePoints(newPoints);

    // Copy new points
    storedPoints() = newPoints;
}


template<class Face>
void Foam::MeshedSurface<Face>::scalePoints(const scalar scaleFactor)
{
    // avoid bad scaling
    if (scaleFactor > 0 && scaleFactor != 1.0)
    {
        pointField newPoints(scaleFactor*this->points());

        // Adapt for new point position
        ParentType::movePoints(newPoints);

        storedPoints() = newPoints;
    }
}


template<class Face>
void Foam::MeshedSurface<Face>::reset
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face> >& faceLst,
    const Xfer<surfZoneList>& zoneLst
)
{
    ParentType::clearOut();

    // Take over new primitive data.
    // Optimized to avoid overwriting data at all
    if (&pointLst)
    {
        storedPoints().transfer(pointLst());
    }

    if (&faceLst)
    {
        storedFaces().transfer(faceLst());
    }

    if (&zoneLst)
    {
        storedZones().transfer(zoneLst());
    }
}


template<class Face>
void Foam::MeshedSurface<Face>::reset
(
    const Xfer<List<point> >& pointLst,
    const Xfer<List<Face> >& faceLst,
    const Xfer<surfZoneList>& zoneLst
)
{
    ParentType::clearOut();

    // Take over new primitive data.
    // Optimized to avoid overwriting data at all
    if (&pointLst)
    {
        storedPoints().transfer(pointLst());
    }

    if (&faceLst)
    {
        storedFaces().transfer(faceLst());
    }

    if (&zoneLst)
    {
        storedZones().transfer(zoneLst());
    }
}


// Remove badly degenerate faces, double faces.
template<class Face>
void Foam::MeshedSurface<Face>::cleanup(const bool verbose)
{
    // merge points (already done for STL, TRI)
    stitchFaces(SMALL, verbose);

    checkFaces(verbose);
    this->checkTopology(verbose);
}


template<class Face>
bool Foam::MeshedSurface<Face>::stitchFaces
(
    const scalar tol,
    const bool verbose
)
{
    pointField& pointLst = this->storedPoints();

    // Merge points
    labelList  pointMap(pointLst.size());
    pointField newPoints(pointLst.size());

    bool hasMerged = mergePoints(pointLst, tol, verbose, pointMap, newPoints);

    if (!hasMerged)
    {
        return false;
    }

    if (verbose)
    {
        Info<< "MeshedSurface::stitchFaces : Renumbering all faces"
            << endl;
    }

    // Set the coordinates to the merged ones
    pointLst.transfer(newPoints);

    List<Face>& faceLst = this->storedFaces();

    List<label> faceMap(faceLst.size());

    // Reset the point labels to the unique points array
    label newFaceI = 0;
    forAll(faceLst, faceI)
    {
        Face& f = faceLst[faceI];
        forAll(f, fp)
        {
            f[fp] = pointMap[f[fp]];
        }

        // for extra safety: collapse face as well
        if (f.collapse() >= 3)
        {
            if (newFaceI != faceI)
            {
                faceLst[newFaceI] = f;
            }
            faceMap[newFaceI] = faceI;
            newFaceI++;
        }
        else if (verbose)
        {
            Pout<< "MeshedSurface::stitchFaces : "
                << "Removing collapsed face " << faceI << endl
                << "    vertices   :" << f << endl;
        }
    }
    pointMap.clear();

    if (newFaceI != faceLst.size())
    {
        if (verbose)
        {
            Pout<< "MeshedSurface::stitchFaces : "
                << "Removed " << faceLst.size() - newFaceI
                << " faces" << endl;
        }
        faceLst.setSize(newFaceI);
        faceMap.setSize(newFaceI);
        remapFaces(faceMap);
    }
    faceMap.clear();

    // Merging points might have changed geometric factors
    ParentType::clearOut();
    return true;
}


// Remove badly degenerate faces and double faces.
template<class Face>
bool Foam::MeshedSurface<Face>::checkFaces
(
    const bool verbose
)
{
    bool changed = false;
    List<Face>& faceLst = this->storedFaces();

    List<label> faceMap(faceLst.size());

    label newFaceI = 0;
    // Detect badly labelled faces and mark degenerate faces
    const label maxPointI = this->points().size() - 1;
    forAll(faceLst, faceI)
    {
        Face& f = faceLst[faceI];

        // avoid degenerate faces
        if (f.collapse() >= 3)
        {
            forAll(f, fp)
            {
                if (f[fp] < 0 || f[fp] > maxPointI)
                {
                    FatalErrorIn("MeshedSurface::checkFaces(bool)")
                        << "face " << f
                        << " uses point indices outside point range 0.."
                    << maxPointI
                        << exit(FatalError);
                }
            }

            faceMap[faceI] = faceI;
            newFaceI++;
        }
        else
        {
            // mark as bad face
            faceMap[faceI] = -1;

            changed = true;
            if (verbose)
            {
                WarningIn
                (
                    "MeshedSurface::checkFaces(bool verbose)"
                )   << "face[" << faceI << "] = " << f
                    << " does not have three unique vertices" << endl;
            }
        }
    }

    // Detect doubled faces
    // do not touch the faces
    const labelListList& fFaces = this->faceFaces();
    newFaceI = 0;
    forAll(faceLst, faceI)
    {
        // skip already collapsed faces:
        if (faceMap[faceI] < 0)
        {
            continue;
        }

        const Face& f = faceLst[faceI];

        // duplicate face check
        bool okay = true;
        const labelList& neighbours = fFaces[faceI];

        // Check if faceNeighbours use same points as this face.
        // Note: discards normal information - sides of baffle are merged.
        forAll(neighbours, neighI)
        {
            const label neiFaceI = neighbours[neighI];

            if (neiFaceI <= faceI || faceMap[neiFaceI] < 0)
            {
                // lower numbered faces already checked
                // skip neighbours that are themselves collapsed
                continue;
            }

            const Face& nei = faceLst[neiFaceI];

            if (f == nei)
            {
                okay = false;

                if (verbose)
                {
                    WarningIn
                    (
                        "MeshedSurface::checkFaces(bool verbose)"
                    )   << "faces share the same vertices:" << nl
                        << "    face[" << faceI << "] : " << f << nl
                        << "    face[" << neiFaceI << "] : " << nei << endl;
                    // printFace(Warning, "    ", f, points());
                    // printFace(Warning, "    ", nei, points());
                }

                break;
            }
        }

        if (okay)
        {
            faceMap[faceI] = faceI;
            newFaceI++;
        }
        else
        {
            faceMap[faceI] = -1;
        }
    }

    // Phase 1: pack
    // Done to keep numbering constant in phase 1

    if (changed || newFaceI < faceLst.size())
    {
        changed = true;

        if (verbose)
        {
            WarningIn
            (
                "MeshedSurface::checkFaces(bool verbose)"
            )   << "Removed " << faceLst.size() - newFaceI
                << " illegal faces." << endl;
        }

        // compress the face list
        newFaceI = 0;
        forAll(faceLst, faceI)
        {
            if (faceMap[faceI] >= 0)
            {
                if (newFaceI != faceI)
                {
                    faceLst[newFaceI] = faceLst[faceI];
                }
                faceMap[newFaceI] = faceI;
                newFaceI++;
            }
        }

        faceLst.setSize(newFaceI);
        remapFaces(faceMap);
    }
    faceMap.clear();

    // Topology can change because of renumbering
    ParentType::clearOut();
    return changed;
}


template<class Face>
Foam::label Foam::MeshedSurface<Face>::triangulate()
{
    return triangulate
    (
        const_cast<List<label>&>(List<label>::null())
    );
}


template<class Face>
Foam::label Foam::MeshedSurface<Face>::triangulate
(
    List<label>& faceMapOut
)
{
    label nTri = 0;
    label maxTri = 0;  // the maximum number of triangles for any single face
    List<Face>& faceLst = this->storedFaces();

    // determine how many triangles will be needed
    forAll(faceLst, faceI)
    {
        const label n = faceLst[faceI].nTriangles();
        if (maxTri < n)
        {
            maxTri = n;
        }
        nTri += n;
    }

    // nothing to do
    if (nTri <= faceLst.size())
    {
        if (&faceMapOut)
        {
            faceMapOut.clear();
        }
        return 0;
    }

    List<Face>  newFaces(nTri);
    List<label> faceMap;

    // reuse storage from optional faceMap
    if (&faceMapOut)
    {
        faceMap.transfer(faceMapOut);
    }
    faceMap.setSize(nTri);

    // remember the number of *additional* faces
    nTri -= faceLst.size();

    if (this->points().empty())
    {
        // triangulate without points
        // simple face triangulation around f[0]
        label newFaceI = 0;
        forAll(faceLst, faceI)
        {
            const Face& f = faceLst[faceI];

            for (label fp = 1; fp < f.size() - 1; ++fp)
            {
                label fp1 = f.fcIndex(fp);

                newFaces[newFaceI] = triFace(f[0], f[fp], f[fp1]);
                faceMap[newFaceI] = faceI;
                newFaceI++;
            }
        }
    }
    else
    {
        // triangulate with points
        List<face> tmpTri(maxTri);

        label newFaceI = 0;
        forAll(faceLst, faceI)
        {
            // 'face' not '<Face>'
            const face& f = faceLst[faceI];

            label nTmp = 0;
            f.triangles(this->points(), nTmp, tmpTri);
            for (label triI = 0; triI < nTmp; triI++)
            {
                newFaces[newFaceI] = Face
                (
                    static_cast<labelUList&>(tmpTri[triI])
                );
                faceMap[newFaceI] = faceI;
                newFaceI++;
            }
        }
    }

    faceLst.transfer(newFaces);
    remapFaces(faceMap);

    // optionally return the faceMap
    if (&faceMapOut)
    {
        faceMapOut.transfer(faceMap);
    }
    faceMap.clear();

    // Topology can change because of renumbering
    ParentType::clearOut();
    return nTri;
}




template<class Face>
Foam::MeshedSurface<Face> Foam::MeshedSurface<Face>::subsetMesh
(
    const labelHashSet& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    const pointField& locPoints = this->localPoints();
    const List<Face>& locFaces  = this->localFaces();


    // Fill pointMap, faceMap
    PatchTools::subsetMap(*this, include, pointMap, faceMap);

    // Create compact coordinate list and forward mapping array
    pointField newPoints(pointMap.size());
    labelList oldToNew(locPoints.size());
    forAll(pointMap, pointI)
    {
        newPoints[pointI] = locPoints[pointMap[pointI]];
        oldToNew[pointMap[pointI]] = pointI;
    }

    // create/copy a new zones list, each zone with zero size
    surfZoneList newZones(this->surfZones());
    forAll(newZones, zoneI)
    {
        newZones[zoneI].size() = 0;
    }

    // Renumber face node labels
    List<Face> newFaces(faceMap.size());
    forAll(faceMap, faceI)
    {
        const label origFaceI = faceMap[faceI];
        newFaces[faceI] = Face(locFaces[origFaceI]);

        // Renumber labels for face
        Face& f = newFaces[faceI];
        forAll(f, fp)
        {
            f[fp] = oldToNew[f[fp]];
        }
    }
    oldToNew.clear();

    // recalculate the zones start/size
    label newFaceI = 0;
    label origEndI = 0;

    // adjust zone sizes
    forAll(newZones, zoneI)
    {
        surfZone& zone = newZones[zoneI];

        // adjust zone start
        zone.start() = newFaceI;
        origEndI += zone.size();

        for (label faceI = newFaceI; faceI < faceMap.size(); ++faceI)
        {
            if (faceMap[faceI] < origEndI)
            {
                ++newFaceI;
            }
            else
            {
                break;
            }
        }

        // adjust zone size
        zone.size() = newFaceI - zone.start();
    }


    // construct a sub-surface
    return MeshedSurface
    (
        xferMove(newPoints),
        xferMove(newFaces),
        xferMove(newZones)
    );
}


template<class Face>
Foam::MeshedSurface<Face> Foam::MeshedSurface<Face>::subsetMesh
(
    const labelHashSet& include
) const
{
    labelList pointMap, faceMap;
    return subsetMesh(include, pointMap, faceMap);
}



template<class Face>
void Foam::MeshedSurface<Face>::transfer
(
    MeshedSurface<Face>& surf
)
{
    reset
    (
        xferMove(surf.storedPoints()),
        xferMove(surf.storedFaces()),
        xferMove(surf.storedZones())
    );
}


template<class Face>
void Foam::MeshedSurface<Face>::transfer
(
    UnsortedMeshedSurface<Face>& surf
)
{
    clear();

    labelList faceMap;
    surfZoneList zoneLst = surf.sortedZones(faceMap);

    if (zoneLst.size() <= 1)
    {
        reset
        (
            xferMove(surf.storedPoints()),
            xferMove(surf.storedFaces()),
            Xfer<surfZoneList>()
        );
    }
    else
    {
        List<Face>& oldFaces = surf.storedFaces();
        List<Face> newFaces(faceMap.size());

        forAll(faceMap, faceI)
        {
            newFaces[faceMap[faceI]].transfer(oldFaces[faceI]);
        }

        reset
        (
            xferMove(surf.storedPoints()),
            xferMove(newFaces),
            xferMove(zoneLst)
        );
    }

    faceMap.clear();
    surf.clear();
}


template<class Face>
Foam::Xfer<Foam::MeshedSurface<Face> > Foam::MeshedSurface<Face>::xfer()
{
    return xferMove(*this);
}


// Read from file, determine format from extension
template<class Face>
bool Foam::MeshedSurface<Face>::read(const fileName& name)
{
    word ext = name.ext();
    if (ext == "gz")
    {
        fileName unzipName = name.lessExt();
        return read(unzipName, unzipName.ext());
    }
    else
    {
        return read(name, ext);
    }
}


// Read from file in given format
template<class Face>
bool Foam::MeshedSurface<Face>::read
(
    const fileName& name,
    const word& ext
)
{
    clear();

    // read via selector mechanism
    transfer(New(name, ext)());
    return true;
}


template<class Face>
void Foam::MeshedSurface<Face>::write
(
    const Time& t,
    const word& surfName
) const
{
    MeshedSurfaceProxy<Face>(*this).write(t, surfName);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::operator=(const MeshedSurface& surf)
{
    clear();

    this->storedPoints() = surf.points();
    this->storedFaces()  = surf.faces();
    this->storedZones()  = surf.surfZones();
}


template<class Face>
Foam::MeshedSurface<Face>::operator Foam::MeshedSurfaceProxy<Face>() const
{
    return MeshedSurfaceProxy<Face>
    (
        this->points(),
        this->faces(),
        this->surfZones()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MeshedSurfaceZones.C"
#include "MeshedSurfaceIO.C"
#include "MeshedSurfaceNew.C"

// ************************************************************************* //
