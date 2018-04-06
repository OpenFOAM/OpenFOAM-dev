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

#include "triSurface.H"
#include "demandDrivenData.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"
#include "boundBox.H"
#include "SortableList.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(triSurface, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::triSurface::triSurfInstance(const Time& d)
{
    fileName foamName(d.caseName() + ".ftr");

    // Search back through the time directories list to find the time
    // closest to and lower than current time

    instantList ts = d.times();
    label i;

    for (i=ts.size()-1; i>=0; i--)
    {
        if (ts[i].value() <= d.timeOutputValue())
        {
            break;
        }
    }

    // Noting that the current directory has already been searched
    // for mesh data, start searching from the previously stored time directory

    if (i>=0)
    {
        for (label j=i; j>=0; j--)
        {
            if (isFile(d.path()/ts[j].name()/typeName/foamName))
            {
                if (debug)
                {
                    Pout<< " triSurface::triSurfInstance(const Time& d)"
                        << "reading " << foamName
                        << " from " << ts[j].name()/typeName
                        << endl;
                }

                return ts[j].name();
            }
        }
    }

    if (debug)
    {
        Pout<< " triSurface::triSurfInstance(const Time& d)"
            << "reading " << foamName
            << " from constant/" << endl;
    }
    return d.constant();
}


Foam::List<Foam::labelledTri> Foam::triSurface::convertToTri
(
    const faceList& faces,
    const label defaultRegion
)
{
    List<labelledTri> triFaces(faces.size());

    forAll(triFaces, facei)
    {
        const face& f = faces[facei];

        if (f.size() != 3)
        {
            FatalErrorInFunction
                << "Face at position " << facei
                << " does not have three vertices:" << f
                << abort(FatalError);
        }

        labelledTri& tri = triFaces[facei];

        tri[0] = f[0];
        tri[1] = f[1];
        tri[2] = f[2];
        tri.region() = defaultRegion;
    }

    return triFaces;
}


Foam::List<Foam::labelledTri> Foam::triSurface::convertToTri
(
    const triFaceList& faces,
    const label defaultRegion
)
{
    List<labelledTri> triFaces(faces.size());

    forAll(triFaces, facei)
    {
        const triFace& f = faces[facei];

        labelledTri& tri = triFaces[facei];

        tri[0] = f[0];
        tri[1] = f[1];
        tri[2] = f[2];
        tri.region() = defaultRegion;
    }

    return triFaces;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::triSurface::printTriangle
(
    Ostream& os,
    const string& pre,
    const labelledTri& f,
    const pointField& points
)
{
    os
        << pre.c_str() << "vertex numbers:"
        << f[0] << ' ' << f[1] << ' ' << f[2] << endl
        << pre.c_str() << "vertex coords :"
        << points[f[0]] << ' ' << points[f[1]] << ' ' << points[f[2]]
        << pre.c_str() << "region        :" << f.region() << endl
        << endl;
}


Foam::string Foam::triSurface::getLineNoComment(IFstream& is)
{
    string line;
    do
    {
        is.getLine(line);
    }
    while ((line.empty() || line[0] == '#') && is.good());

    return line;
}


// Remove non-triangles, double triangles.
void Foam::triSurface::checkTriangles(const bool verbose)
{
    // Simple check on indices ok.
    const label maxPointi = points().size() - 1;

    forAll(*this, facei)
    {
        const triSurface::FaceType& f = (*this)[facei];

        forAll(f, fp)
        {
            if (f[fp] < 0 || f[fp] > maxPointi)
            {
                FatalErrorInFunction
                    << "triangle " << f
                    << " uses point indices outside point range 0.."
                    << maxPointi
                    << exit(FatalError);
            }
        }
    }

    // Two phase process
    //   1. mark invalid faces
    //   2. pack
    // Done to keep numbering constant in phase 1

    // List of valid triangles
    boolList valid(size(), true);
    bool hasInvalid = false;

    forAll(*this, facei)
    {
        const labelledTri& f = (*this)[facei];

        if ((f[0] == f[1]) || (f[0] == f[2]) || (f[1] == f[2]))
        {
            // 'degenerate' triangle check
            valid[facei] = false;
            hasInvalid = true;

            if (verbose)
            {
                WarningInFunction
                    << "triangle " << facei
                    << " does not have three unique vertices:\n";
                printTriangle(Warning, "    ", f, points());
            }
        }
        else
        {
            // duplicate triangle check
            const labelList& fEdges = faceEdges()[facei];

            // Check if faceNeighbours use same points as this face.
            // Note: discards normal information - sides of baffle are merged.

            forAll(fEdges, fp)
            {
                const labelList& eFaces = edgeFaces()[fEdges[fp]];

                forAll(eFaces, i)
                {
                    label neighbour = eFaces[i];

                    if (neighbour > facei)
                    {
                        // lower numbered faces already checked
                        const labelledTri& n = (*this)[neighbour];

                        if
                        (
                            ((f[0] == n[0]) || (f[0] == n[1]) || (f[0] == n[2]))
                         && ((f[1] == n[0]) || (f[1] == n[1]) || (f[1] == n[2]))
                         && ((f[2] == n[0]) || (f[2] == n[1]) || (f[2] == n[2]))
                        )
                        {
                            valid[facei] = false;
                            hasInvalid = true;

                            if (verbose)
                            {
                                WarningInFunction
                                    << "triangles share the same vertices:\n"
                                    << "    face 1 :" << facei << endl;
                                printTriangle(Warning, "    ", f, points());

                                Warning
                                    << endl
                                    << "    face 2 :"
                                    << neighbour << endl;
                                printTriangle(Warning, "    ", n, points());
                            }

                            break;
                        }
                    }
                }
            }
        }
    }

    if (hasInvalid)
    {
        // Pack
        label newFacei = 0;
        forAll(*this, facei)
        {
            if (valid[facei])
            {
                const labelledTri& f = (*this)[facei];
                (*this)[newFacei++] = f;
            }
        }

        if (verbose)
        {
            WarningInFunction
                << "Removing " << size() - newFacei
                << " illegal faces." << endl;
        }
        (*this).setSize(newFacei);

        // Topology can change because of renumbering
        clearOut();
    }
}


// Check/fix edges with more than two triangles
void Foam::triSurface::checkEdges(const bool verbose)
{
    const labelListList& eFaces = edgeFaces();

    forAll(eFaces, edgeI)
    {
        const labelList& myFaces = eFaces[edgeI];

        if (myFaces.empty())
        {
            FatalErrorInFunction
                << "Edge " << edgeI << " with vertices " << edges()[edgeI]
                << " has no edgeFaces"
                << exit(FatalError);
        }
        else if (myFaces.size() > 2 && verbose)
        {
            WarningInFunction
                << "Edge " << edgeI << " with vertices " << edges()[edgeI]
                << " has more than 2 faces connected to it : " << myFaces
                << endl;
        }
    }
}


// Read triangles, points from Istream
bool Foam::triSurface::read(Istream& is)
{
    is  >> patches_ >> storedPoints() >> storedFaces();

    return true;
}


// Read from file in given format
bool Foam::triSurface::read
(
    const fileName& name,
    const word& ext,
    const bool check
)
{
    if (check && !exists(name))
    {
        FatalErrorInFunction
            << "Cannot read " << name << exit(FatalError);
    }

    if (ext == "gz")
    {
        fileName unzipName = name.lessExt();

        // Do not check for existence. Let IFstream do the unzipping.
        return read(unzipName, unzipName.ext(), false);
    }
    else if (ext == "ftr")
    {
        return read(IFstream(name)());
    }
    else if (ext == "stl")
    {
        return readSTL(name);
    }
    else if (ext == "stlb")
    {
        return readSTLBINARY(name);
    }
    else if (ext == "gts")
    {
        return readGTS(name);
    }
    else if (ext == "obj")
    {
        return readOBJ(name);
    }
    else if (ext == "off")
    {
        return readOFF(name);
    }
    else if (ext == "tri")
    {
        return readTRI(name);
    }
    else if (ext == "ac")
    {
        return readAC(name);
    }
    else if (ext == "nas")
    {
        return readNAS(name);
    }
    else if (ext == "vtk")
    {
        return readVTK(name);
    }
    else
    {
        FatalErrorInFunction
            << "unknown file extension " << ext
            << ". Supported extensions are '.ftr', '.stl', '.stlb', '.gts'"
            << ", '.obj', '.ac', '.off', '.nas', '.tri' and '.vtk'"
            << exit(FatalError);

        return false;
    }
}


// Write to file in given format
void Foam::triSurface::write
(
    const fileName& name,
    const word& ext,
    const bool sort
) const
{
    if (ext == "ftr")
    {
        return write(OFstream(name)());
    }
    else if (ext == "stl")
    {
        return writeSTLASCII(sort, OFstream(name)());
    }
    else if (ext == "stlb")
    {
        ofstream outFile(name.c_str(), std::ios::binary);

        writeSTLBINARY(outFile);
    }
    else if (ext == "gts")
    {
        return writeGTS(sort, OFstream(name)());
    }
    else if (ext == "obj")
    {
        writeOBJ(sort, OFstream(name)());
    }
    else if (ext == "off")
    {
        writeOFF(sort, OFstream(name)());
    }
    else if (ext == "vtk")
    {
        writeVTK(sort, OFstream(name)());
    }
    else if (ext == "tri")
    {
        writeTRI(sort, OFstream(name)());
    }
    else if (ext == "ac")
    {
        writeAC(OFstream(name)());
    }
    else if (ext == "smesh")
    {
        writeSMESH(sort, OFstream(name)());
    }
    else
    {
        FatalErrorInFunction
            << "unknown file extension " << ext
            << " for file " << name
            << ". Supported extensions are '.ftr', '.stl', '.stlb', "
            << "'.gts', '.obj', '.vtk'"
            << ", '.off', '.dx', '.smesh', '.ac' and '.tri'"
            << exit(FatalError);
    }
}


// Returns patch info. Sets faceMap to the indexing according to patch
// numbers. Patch numbers start at 0.
Foam::surfacePatchList Foam::triSurface::calcPatches(labelList& faceMap) const
{
    // Sort according to region numbers of labelledTri
    SortableList<label> sortedRegion(size());

    forAll(sortedRegion, facei)
    {
        sortedRegion[facei] = operator[](facei).region();
    }
    sortedRegion.sort();

    faceMap = sortedRegion.indices();

    // Extend regions
    label maxRegion = patches_.size()-1;    // for non-compacted regions

    if (faceMap.size())
    {
        maxRegion = max
        (
            maxRegion,
            operator[](faceMap.last()).region()
        );
    }

    // Get new region list
    surfacePatchList newPatches(maxRegion + 1);

    // Fill patch sizes
    forAll(*this, facei)
    {
        label region = operator[](facei).region();

        newPatches[region].size()++;
    }


    // Fill rest of patch info

    label startFacei = 0;
    forAll(newPatches, newPatchi)
    {
        surfacePatch& newPatch = newPatches[newPatchi];

        newPatch.index() = newPatchi;

        label oldPatchi = newPatchi;

        // start of patch
        newPatch.start() = startFacei;


        // Take over any information from existing patches
        if ((oldPatchi < patches_.size()) && (patches_[oldPatchi].name() != ""))
        {
            newPatch.name() = patches_[oldPatchi].name();
        }
        else
        {
            newPatch.name() = word("patch") + name(newPatchi);
        }

        if
        (
            (oldPatchi < patches_.size())
         && (patches_[oldPatchi].geometricType() != "")
        )
        {
            newPatch.geometricType() = patches_[oldPatchi].geometricType();
        }
        else
        {
            newPatch.geometricType() = "empty";
        }

        startFacei += newPatch.size();
    }

    return newPatches;
}


void Foam::triSurface::setDefaultPatches()
{
    labelList faceMap;

    // Get names, types and sizes
    surfacePatchList newPatches(calcPatches(faceMap));

    // Take over names and types (but not size)
    patches_.setSize(newPatches.size());

    forAll(newPatches, patchi)
    {
        patches_[patchi].index() = patchi;
        patches_[patchi].name() = newPatches[patchi].name();
        patches_[patchi].geometricType() = newPatches[patchi].geometricType();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurface::triSurface()
:
    ParentType(List<Face>(), pointField()),
    patches_(0),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{}



Foam::triSurface::triSurface
(
    const List<labelledTri>& triangles,
    const geometricSurfacePatchList& patches,
    const pointField& points
)
:
    ParentType(triangles, points),
    patches_(patches),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{}


Foam::triSurface::triSurface
(
    List<labelledTri>& triangles,
    const geometricSurfacePatchList& patches,
    pointField& points,
    const bool reuse
)
:
    ParentType(triangles, points, reuse),
    patches_(patches),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{}


Foam::triSurface::triSurface
(
    const Xfer<List<labelledTri>>& triangles,
    const geometricSurfacePatchList& patches,
    const Xfer<List<point>>& points
)
:
    ParentType(triangles, points),
    patches_(patches),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{}


Foam::triSurface::triSurface
(
    const List<labelledTri>& triangles,
    const pointField& points
)
:
    ParentType(triangles, points),
    patches_(),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{
    setDefaultPatches();
}


Foam::triSurface::triSurface
(
    const triFaceList& triangles,
    const pointField& points
)
:
    ParentType(convertToTri(triangles, 0), points),
    patches_(),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{
    setDefaultPatches();
}


Foam::triSurface::triSurface(const fileName& name)
:
    ParentType(List<Face>(), pointField()),
    patches_(),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{
    word ext = name.ext();

    read(name, ext);

    setDefaultPatches();
}


Foam::triSurface::triSurface(Istream& is)
:
    ParentType(List<Face>(), pointField()),
    patches_(),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{
    read(is);

    setDefaultPatches();
}


Foam::triSurface::triSurface(const Time& d)
:
    ParentType(List<Face>(), pointField()),
    patches_(),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{
    fileName foamFile(d.caseName() + ".ftr");

    fileName foamPath(d.path()/triSurfInstance(d)/typeName/foamFile);

    IFstream foamStream(foamPath);

    read(foamStream);

    setDefaultPatches();
}


Foam::triSurface::triSurface(const triSurface& ts)
:
    ParentType(ts, ts.points()),
    patches_(ts.patches()),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurface::~triSurface()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triSurface::clearTopology()
{
    ParentType::clearTopology();
    deleteDemandDrivenData(sortedEdgeFacesPtr_);
    deleteDemandDrivenData(edgeOwnerPtr_);
}


void Foam::triSurface::clearPatchMeshAddr()
{
    ParentType::clearPatchMeshAddr();
}


void Foam::triSurface::clearOut()
{
    ParentType::clearOut();

    clearTopology();
    clearPatchMeshAddr();
}


const Foam::labelListList& Foam::triSurface::sortedEdgeFaces() const
{
    if (!sortedEdgeFacesPtr_)
    {
        calcSortedEdgeFaces();
    }

    return *sortedEdgeFacesPtr_;
}


const Foam::labelList& Foam::triSurface::edgeOwner() const
{
    if (!edgeOwnerPtr_)
    {
        calcEdgeOwner();
    }

    return *edgeOwnerPtr_;
}


void Foam::triSurface::movePoints(const pointField& newPoints)
{
    // Remove all geometry dependent data
    deleteDemandDrivenData(sortedEdgeFacesPtr_);

    // Adapt for new point position
    ParentType::movePoints(newPoints);

    // Copy new points
    storedPoints() = newPoints;
}


void Foam::triSurface::scalePoints(const scalar scaleFactor)
{
    // avoid bad scaling
    if (scaleFactor > 0 && scaleFactor != 1.0)
    {
        // Remove all geometry dependent data
        clearTopology();

        // Adapt for new point position
        ParentType::movePoints(pointField());

        storedPoints() *= scaleFactor;
    }
}


// Remove non-triangles, double triangles.
void Foam::triSurface::cleanup(const bool verbose)
{
    // Merge points (already done for STL, TRI)
    stitchTriangles(small, verbose);

    // Merging points might have changed geometric factors
    clearOut();

    checkTriangles(verbose);

    checkEdges(verbose);
}


// Finds area, starting at facei, delimited by borderEdge. Marks all visited
// faces (from face-edge-face walk) with currentZone.
void Foam::triSurface::markZone
(
    const boolList& borderEdge,
    const label facei,
    const label currentZone,
    labelList& faceZone
) const
{
    // List of faces whose faceZone has been set.
    labelList changedFaces(1, facei);

    while (true)
    {
        // Pick up neighbours of changedFaces
        DynamicList<label> newChangedFaces(2*changedFaces.size());

        forAll(changedFaces, i)
        {
            label facei = changedFaces[i];

            const labelList& fEdges = faceEdges()[facei];

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];

                if (!borderEdge[edgeI])
                {
                    const labelList& eFaces = edgeFaces()[edgeI];

                    forAll(eFaces, j)
                    {
                        label nbrFacei = eFaces[j];

                        if (faceZone[nbrFacei] == -1)
                        {
                            faceZone[nbrFacei] = currentZone;
                            newChangedFaces.append(nbrFacei);
                        }
                        else if (faceZone[nbrFacei] != currentZone)
                        {
                            FatalErrorInFunction
                                << "Zones " << faceZone[nbrFacei]
                                << " at face " << nbrFacei
                                << " connects to zone " << currentZone
                                << " at face " << facei
                                << abort(FatalError);
                        }
                    }
                }
            }
        }

        if (newChangedFaces.empty())
        {
            break;
        }

        changedFaces.transfer(newChangedFaces);
    }
}


// Finds areas delimited by borderEdge (or 'real' edges).
// Fills faceZone accordingly
Foam::label Foam::triSurface::markZones
(
    const boolList& borderEdge,
    labelList& faceZone
) const
{
    faceZone.setSize(size());
    faceZone = -1;

    if (borderEdge.size() != nEdges())
    {
        FatalErrorInFunction
            << "borderEdge boolList not same size as number of edges" << endl
            << "borderEdge:" << borderEdge.size() << endl
            << "nEdges    :" << nEdges()
            << exit(FatalError);
    }

    label zoneI = 0;

    for (label startFacei = 0;; zoneI++)
    {
        // Find first non-coloured face
        for (; startFacei < size(); startFacei++)
        {
            if (faceZone[startFacei] == -1)
            {
                break;
            }
        }

        if (startFacei >= size())
        {
            break;
        }

        faceZone[startFacei] = zoneI;

        markZone(borderEdge, startFacei, zoneI, faceZone);
    }

    return zoneI;
}


void Foam::triSurface::subsetMeshMap
(
    const boolList& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    const List<labelledTri>& locFaces = localFaces();

    label facei = 0;
    label pointi = 0;

    faceMap.setSize(locFaces.size());

    pointMap.setSize(nPoints());

    boolList pointHad(nPoints(), false);

    forAll(include, oldFacei)
    {
        if (include[oldFacei])
        {
            // Store new faces compact
            faceMap[facei++] = oldFacei;

            // Renumber labels for face
            const triSurface::FaceType& f = locFaces[oldFacei];

            forAll(f, fp)
            {
                label labI = f[fp];
                if (!pointHad[labI])
                {
                    pointHad[labI] = true;
                    pointMap[pointi++] = labI;
                }
            }
        }
    }

    // Trim
    faceMap.setSize(facei);
    pointMap.setSize(pointi);
}


Foam::triSurface Foam::triSurface::subsetMesh
(
    const boolList& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    const pointField& locPoints = localPoints();
    const List<labelledTri>& locFaces = localFaces();

    // Fill pointMap, faceMap
    subsetMeshMap(include, pointMap, faceMap);


    // Create compact coordinate list and forward mapping array
    pointField newPoints(pointMap.size());
    labelList oldToNew(locPoints.size());
    forAll(pointMap, pointi)
    {
        newPoints[pointi] = locPoints[pointMap[pointi]];
        oldToNew[pointMap[pointi]] = pointi;
    }

    // Renumber triangle node labels and compact
    List<labelledTri> newTriangles(faceMap.size());

    forAll(faceMap, facei)
    {
        // Get old vertex labels
        const labelledTri& tri = locFaces[faceMap[facei]];

        newTriangles[facei][0] = oldToNew[tri[0]];
        newTriangles[facei][1] = oldToNew[tri[1]];
        newTriangles[facei][2] = oldToNew[tri[2]];
        newTriangles[facei].region() = tri.region();
    }

    // Construct subsurface
    return triSurface(newTriangles, patches(), newPoints, true);
}


Foam::faceList Foam::triSurface::faces() const
{
    faceList faces(size());

    forAll(*this, facei)
    {
        faces[facei] = operator[](facei).triFaceFace();
    }

    return faces;
}


Foam::scalar Foam::triSurface::pointNormalWeight
(
    const triFace& f,
    const label pi,
    const vector& fa,
    const pointField& points
) const
{
    const label index = findIndex(f, pi);

    if (index == -1)
    {
        FatalErrorInFunction
            << "Point not in face" << abort(FatalError);
    }

    const vector e1 = points[f[index]] - points[f[f.fcIndex(index)]];
    const vector e2 = points[f[index]] - points[f[f.rcIndex(index)]];

    return mag(fa)/(magSqr(e1)*magSqr(e2) + vSmall);
}


Foam::tmp<Foam::vectorField> Foam::triSurface::pointNormals2() const
{
    // Weighted average of normals of faces attached to the vertex
    // Weight = fA / (mag(e0)^2 * mag(e1)^2);

    tmp<vectorField> tpointNormals
    (
        new vectorField(nPoints(), Zero)
    );
    vectorField& pointNormals = tpointNormals.ref();

    const pointField& points = this->points();
    const labelListList& pointFaces = this->pointFaces();
    const labelList& meshPoints = this->meshPoints();

    forAll(pointFaces, pi)
    {
        const labelList& pFaces = pointFaces[pi];

        forAll(pFaces, fi)
        {
            const label facei = pFaces[fi];
            const triFace& f = operator[](facei);

            const vector fa = f.area(points);
            const scalar weight =
                pointNormalWeight(f, meshPoints[pi], fa, points);

            pointNormals[pi] += weight*fa;
        }

        pointNormals[pi] /= mag(pointNormals[pi]) + vSmall;
    }

    return tpointNormals;
}


void Foam::triSurface::write
(
    const fileName& name,
    const bool sortByRegion
) const
{
    write(name, name.ext(), sortByRegion);
}


void Foam::triSurface::write(Ostream& os) const
{
    os  << patches() << endl;

    //Note: Write with global point numbering
    os  << points() << nl
        << static_cast<const List<labelledTri>&>(*this) << endl;

    // Check state of Ostream
    os.check("triSurface::write(Ostream&)");
}


void Foam::triSurface::write(const Time& d) const
{
    fileName foamFile(d.caseName() + ".ftr");

    fileName foamPath(d.path()/triSurfInstance(d)/typeName/foamFile);

    OFstream foamStream(foamPath);

    write(foamStream);
}


void Foam::triSurface::writeStats(Ostream& os) const
{
    // Unfortunately nPoints constructs meshPoints() so do compact version
    // ourselves.
    PackedBoolList pointIsUsed(points().size());

    label nPoints = 0;
    boundBox bb = boundBox::invertedBox;

    forAll(*this, facei)
    {
        const triSurface::FaceType& f = operator[](facei);

        forAll(f, fp)
        {
            label pointi = f[fp];
            if (pointIsUsed.set(pointi, 1))
            {
                bb.min() = ::Foam::min(bb.min(), points()[pointi]);
                bb.max() = ::Foam::max(bb.max(), points()[pointi]);
                nPoints++;
            }
        }
    }

    os  << "Triangles    : " << size() << endl
        << "Vertices     : " << nPoints << endl
        << "Bounding Box : " << bb << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::triSurface::operator=(const triSurface& ts)
{
    List<labelledTri>::operator=(ts);
    clearOut();
    storedPoints() = ts.points();
    patches_ = ts.patches();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const triSurface& sm)
{
    sm.write(os);
    return os;
}


// ************************************************************************* //
