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

#include "triSurfaceMesh.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"
#include "EdgeMap.H"
#include "triSurfaceFields.H"
#include "Time.H"
#include "PatchTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(triSurfaceMesh, 0);
    addToRunTimeSelectionTable(searchableSurface, triSurfaceMesh, dict);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::triSurfaceMesh::checkFile
(
    const regIOobject& io,
    const bool isGlobal
)
{
    const fileName fName
    (
        isGlobal
      ? io.globalFilePath(typeName)
      : io.localFilePath(typeName)
    );
    if (fName.empty())
    {
        FatalErrorInFunction
            << "Cannot find triSurfaceMesh starting from "
            << io.objectPath() << exit(FatalError);
    }

    return fName;
}


Foam::fileName Foam::triSurfaceMesh::relativeFilePath
(
    const regIOobject& io,
    const fileName& f,
    const bool isGlobal
)
{
    fileName fName(f);
    fName.expand();
    if (!fName.isAbsolute())
    {
        // Is the specified file:
        // - local to the cwd?
        // - local to the case dir?
        // - or just another name?
        fName = fileHandler().filePath
        (
            isGlobal,
            IOobject(io, fName),
            word::null
        );
    }
    return fName;
}


Foam::fileName Foam::triSurfaceMesh::checkFile
(
    const regIOobject& io,
    const dictionary& dict,
    const bool isGlobal
)
{
    fileName dictFName, fName;

    if (dict.readIfPresent("file", dictFName, false, false))
    {
        fName = relativeFilePath(io, dictFName, isGlobal);

        if (!exists(fName))
        {
            FatalErrorInFunction
                << "Cannot find triSurfaceMesh at " << io.path(dictFName)
                << exit(FatalError);
        }
    }
    else
    {
        fName =
        (
            isGlobal
          ? io.globalFilePath(typeName)
          : io.localFilePath(typeName)
        );

        if (!exists(fName))
        {
            FatalErrorInFunction
                << "Cannot find triSurfaceMesh starting from "
                << io.objectPath() << exit(FatalError);
        }
    }

    return fName;
}


bool Foam::triSurfaceMesh::addFaceToEdge
(
    const edge& e,
    EdgeMap<label>& facesPerEdge
)
{
    EdgeMap<label>::iterator eFnd = facesPerEdge.find(e);
    if (eFnd != facesPerEdge.end())
    {
        if (eFnd() == 2)
        {
            return false;
        }
        eFnd()++;
    }
    else
    {
        facesPerEdge.insert(e, 1);
    }
    return true;
}


bool Foam::triSurfaceMesh::isSurfaceClosed() const
{
    const pointField& pts = triSurface::points();

    // Construct pointFaces. Let's hope surface has compact point
    // numbering ...
    labelListList pointFaces;
    invertManyToMany(pts.size(), *this, pointFaces);

    // Loop over all faces surrounding point. Count edges emanating from point.
    // Every edge should be used by two faces exactly.
    // To prevent doing work twice per edge only look at edges to higher
    // point
    EdgeMap<label> facesPerEdge(100);
    forAll(pointFaces, pointi)
    {
        const labelList& pFaces = pointFaces[pointi];

        facesPerEdge.clear();
        forAll(pFaces, i)
        {
            const triSurface::FaceType& f = triSurface::operator[](pFaces[i]);
            label fp = findIndex(f, pointi);

            // Something weird: if I expand the code of addFaceToEdge in both
            // below instances it gives a segmentation violation on some
            // surfaces. Compiler (4.3.2) problem?


            // Forward edge
            label nextPointi = f[f.fcIndex(fp)];

            if (nextPointi > pointi)
            {
                bool okFace = addFaceToEdge
                (
                    edge(pointi, nextPointi),
                    facesPerEdge
                );

                if (!okFace)
                {
                    return false;
                }
            }
            // Reverse edge
            label prevPointi = f[f.rcIndex(fp)];

            if (prevPointi > pointi)
            {
                bool okFace = addFaceToEdge
                (
                    edge(pointi, prevPointi),
                    facesPerEdge
                );

                if (!okFace)
                {
                    return false;
                }
            }
        }

        // Check for any edges used only once.
        forAllConstIter(EdgeMap<label>, facesPerEdge, iter)
        {
            if (iter() != 2)
            {
                return false;
            }
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io, const triSurface& s)
:
    searchableSurface(io),
    objectRegistry
    (
        IOobject
        (
            io.name(),
            io.instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface(s),
    triSurfaceRegionSearch(s),
    minQuality_(-1),
    surfaceClosed_(-1)
{
    const pointField& pts = triSurface::points();

    bounds() = boundBox(pts);
}


Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io)
:
    // Find instance for triSurfaceMesh
    searchableSurface(io),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            searchableSurface::instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface(checkFile(static_cast<const searchableSurface&>(*this), true)),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this)),
    minQuality_(-1),
    surfaceClosed_(-1)
{
    const pointField& pts = triSurface::points();

    bounds() = boundBox(pts);
}


Foam::triSurfaceMesh::triSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            searchableSurface::instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface
    (
        checkFile(static_cast<const searchableSurface&>(*this), dict, true)
    ),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this), dict),
    minQuality_(-1),
    surfaceClosed_(-1)
{
    // Reading from supplied file name instead of objectPath/filePath
    if (dict.readIfPresent("file", fName_, false, false))
    {
        fName_ = relativeFilePath
        (
            static_cast<const searchableSurface&>(*this),
            fName_,
            true
        );
    }

    scalar scaleFactor = 0;

    // Allow rescaling of the surface points
    // eg, CAD geometries are often done in millimeters
    if (dict.readIfPresent("scale", scaleFactor) && scaleFactor > 0)
    {
        Info<< searchableSurface::name() << " : using scale " << scaleFactor
            << endl;
        triSurface::scalePoints(scaleFactor);
    }

    const pointField& pts = triSurface::points();

    bounds() = boundBox(pts);

    // Have optional minimum quality for normal calculation
    if (dict.readIfPresent("minQuality", minQuality_) && minQuality_ > 0)
    {
        Info<< searchableSurface::name()
            << " : ignoring triangles with quality < "
            << minQuality_ << " for normals calculation." << endl;
    }
}


Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io, const bool isGlobal)
:
    // Find instance for triSurfaceMesh
    searchableSurface(io),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            searchableSurface::instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface
    (
        checkFile(static_cast<const searchableSurface&>(*this), isGlobal)
    ),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this)),
    minQuality_(-1),
    surfaceClosed_(-1)
{
    const pointField& pts = triSurface::points();

    bounds() = boundBox(pts);
}


Foam::triSurfaceMesh::triSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict,
    const bool isGlobal
)
:
    searchableSurface(io),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            searchableSurface::instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface
    (
        checkFile(static_cast<const searchableSurface&>(*this), dict, isGlobal)
    ),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this), dict),
    minQuality_(-1),
    surfaceClosed_(-1)
{
    // Reading from supplied file name instead of objectPath/filePath
    if (dict.readIfPresent("file", fName_, false, false))
    {
        fName_ = relativeFilePath
        (
            static_cast<const searchableSurface&>(*this),
            fName_,
            isGlobal
        );
    }

    scalar scaleFactor = 0;

    // Allow rescaling of the surface points
    // eg, CAD geometries are often done in millimeters
    if (dict.readIfPresent("scale", scaleFactor) && scaleFactor > 0)
    {
        Info<< searchableSurface::name() << " : using scale " << scaleFactor
            << endl;
        triSurface::scalePoints(scaleFactor);
    }

    const pointField& pts = triSurface::points();

    bounds() = boundBox(pts);

    // Have optional minimum quality for normal calculation
    if (dict.readIfPresent("minQuality", minQuality_) && minQuality_ > 0)
    {
        Info<< searchableSurface::name()
            << " : ignoring triangles with quality < "
            << minQuality_ << " for normals calculation." << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceMesh::~triSurfaceMesh()
{
    clearOut();
}


void Foam::triSurfaceMesh::clearOut()
{
    triSurfaceRegionSearch::clearOut();
    edgeTree_.clear();
    triSurface::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::triSurfaceMesh::coordinates() const
{
    tmp<pointField> tPts(new pointField(8));
    pointField& pt = tPts.ref();

    // Use copy to calculate face centres so they don't get stored
    pt = PrimitivePatch<triSurface::FaceType, SubList, const pointField&>
    (
        SubList<triSurface::FaceType>(*this, triSurface::size()),
        triSurface::points()
    ).faceCentres();

    return tPts;
}


void Foam::triSurfaceMesh::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres = coordinates();
    radiusSqr.setSize(size());
    radiusSqr = 0.0;

    const pointField& pts = triSurface::points();

    forAll(*this, facei)
    {
        const labelledTri& f = triSurface::operator[](facei);
        const point& fc = centres[facei];
        forAll(f, fp)
        {
            const point& pt = pts[f[fp]];
            radiusSqr[facei] = max(radiusSqr[facei], Foam::magSqr(fc-pt));
        }
    }

    // Add a bit to make sure all points are tested inside
    radiusSqr += Foam::sqr(small);
}


Foam::tmp<Foam::pointField> Foam::triSurfaceMesh::points() const
{
    return triSurface::points();
}


bool Foam::triSurfaceMesh::overlaps(const boundBox& bb) const
{
     const indexedOctree<treeDataTriSurface>& octree = tree();

     labelList indices = octree.findBox(treeBoundBox(bb));

     return !indices.empty();
}


void Foam::triSurfaceMesh::movePoints(const pointField& newPoints)
{
    triSurfaceRegionSearch::clearOut();
    edgeTree_.clear();
    triSurface::movePoints(newPoints);
}


const Foam::indexedOctree<Foam::treeDataEdge>&
Foam::triSurfaceMesh::edgeTree() const
{
    if (edgeTree_.empty())
    {
        // Boundary edges
        labelList bEdges
        (
            identity
            (
                nEdges()
               -nInternalEdges()
            )
          + nInternalEdges()
        );

        treeBoundBox bb(Zero, Zero);

        if (bEdges.size())
        {
            label nPoints;
            PatchTools::calcBounds
            (
                *this,
                bb,
                nPoints
            );

            // Random number generator. Bit dodgy since not exactly random ;-)
            Random rndGen(65431);

            // Slightly extended bb. Slightly off-centred just so on symmetric
            // geometry there are less face/edge aligned items.

            bb = bb.extend(rndGen, 1e-4);
            bb.min() -= point(rootVSmall, rootVSmall, rootVSmall);
            bb.max() += point(rootVSmall, rootVSmall, rootVSmall);
        }

        scalar oldTol = indexedOctree<treeDataEdge>::perturbTol();
        indexedOctree<treeDataEdge>::perturbTol() = tolerance();

        edgeTree_.reset
        (
            new indexedOctree<treeDataEdge>
            (
                treeDataEdge
                (
                    false,          // cachebb
                    edges(),        // edges
                    localPoints(),  // points
                    bEdges          // selected edges
                ),
                bb,                 // bb
                maxTreeDepth(),     // maxLevel
                10,                 // leafsize
                3.0                 // duplicity
            )
        );

        indexedOctree<treeDataEdge>::perturbTol() = oldTol;
    }
    return edgeTree_();
}


const Foam::wordList& Foam::triSurfaceMesh::regions() const
{
    if (regions_.empty())
    {
        regions_.setSize(patches().size());
        forAll(regions_, regionI)
        {
            regions_[regionI] = patches()[regionI].name();
        }
    }
    return regions_;
}


bool Foam::triSurfaceMesh::hasVolumeType() const
{
    if (surfaceClosed_ == -1)
    {
        if (isSurfaceClosed())
        {
            surfaceClosed_ = 1;
        }
        else
        {
            surfaceClosed_ = 0;
        }
    }

    return surfaceClosed_ == 1;
}


void Foam::triSurfaceMesh::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    triSurfaceSearch::findNearest(samples, nearestDistSqr, info);
}


void Foam::triSurfaceMesh::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    const labelList& regionIndices,
    List<pointIndexHit>& info
) const
{
    triSurfaceRegionSearch::findNearest
    (
        samples,
        nearestDistSqr,
        regionIndices,
        info
    );
}


void Foam::triSurfaceMesh::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    triSurfaceSearch::findLine(start, end, info);
}


void Foam::triSurfaceMesh::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    triSurfaceSearch::findLineAny(start, end, info);
}


void Foam::triSurfaceMesh::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    triSurfaceSearch::findLineAll(start, end, info);
}


void Foam::triSurfaceMesh::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.setSize(info.size());
    forAll(info, i)
    {
        if (info[i].hit())
        {
            region[i] = triSurface::operator[](info[i].index()).region();
        }
        else
        {
            region[i] = -1;
        }
    }
}


void Foam::triSurfaceMesh::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    const triSurface& s = *this;
    const pointField& pts = s.points();

    normal.setSize(info.size());

    if (minQuality_ >= 0)
    {
        // Make sure we don't use triangles with low quality since
        // normal is not reliable.

        const labelListList& faceFaces = s.faceFaces();

        forAll(info, i)
        {
            if (info[i].hit())
            {
                label facei = info[i].index();
                normal[i] = s[facei].area(pts);

                scalar qual = s[facei].tri(pts).quality();

                if (qual < minQuality_)
                {
                    // Search neighbouring triangles
                    const labelList& fFaces = faceFaces[facei];

                    forAll(fFaces, j)
                    {
                        label nbrI = fFaces[j];
                        scalar nbrQual = s[nbrI].tri(pts).quality();
                        if (nbrQual > qual)
                        {
                            qual = nbrQual;
                            normal[i] = s[nbrI].area(pts);
                        }
                    }
                }

                normal[i] /= mag(normal[i]) + vSmall;
            }
            else
            {
                // Set to what?
                normal[i] = Zero;
            }
        }
    }
    else
    {
        forAll(info, i)
        {
            if (info[i].hit())
            {
                label facei = info[i].index();
                // Cached:
                //normal[i] = faceNormals()[facei];

                // Uncached
                normal[i] = s[facei].normal(pts);
            }
            else
            {
                // Set to what?
                normal[i] = Zero;
            }
        }
    }
}


void Foam::triSurfaceMesh::setField(const labelList& values)
{
    autoPtr<triSurfaceLabelField> fldPtr
    (
        new triSurfaceLabelField
        (
            IOobject
            (
                "values",
                objectRegistry::time().timeName(),  // instance
                "triSurface",                       // local
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimless,
            labelField(values)
        )
    );

    // Store field on triMesh
    fldPtr.ptr()->store();
}


void Foam::triSurfaceMesh::getField
(
    const List<pointIndexHit>& info,
    labelList& values
) const
{
    if (foundObject<triSurfaceLabelField>("values"))
    {
        values.setSize(info.size());

        const triSurfaceLabelField& fld = lookupObject<triSurfaceLabelField>
        (
            "values"
        );

        forAll(info, i)
        {
            if (info[i].hit())
            {
                values[i] = fld[info[i].index()];
            }
        }
    }
}


void Foam::triSurfaceMesh::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    volType.setSize(points.size());

    scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = tolerance();

    forAll(points, pointi)
    {
        const point& pt = points[pointi];

        if (!tree().bb().contains(pt))
        {
            // Have to calculate directly as outside the octree
            volType[pointi] = tree().shapes().getVolumeType(tree(), pt);
        }
        else
        {
            // - use cached volume type per each tree node
            volType[pointi] = tree().getVolumeType(pt);
        }
    }

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


bool Foam::triSurfaceMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    fileName fullPath;
    if (fName_.size())
    {
        // Override file name

        fullPath = fName_;

        fullPath.expand();
        if (!fullPath.isAbsolute())
        {
            // Add directory from regIOobject
            fullPath = searchableSurface::objectPath().path()/fullPath;
        }
    }
    else
    {
        fullPath = searchableSurface::objectPath();
    }

    if (!mkDir(fullPath.path()))
    {
        return false;
    }

    triSurface::write(fullPath);

    if (!isFile(fullPath))
    {
        return false;
    }

    return true;
}


// ************************************************************************* //
