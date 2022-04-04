/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "localPointRegion.H"
#include "syncTools.H"
#include "polyMesh.H"
#include "polyTopoChangeMap.H"
#include "globalIndex.H"
#include "indirectPrimitivePatch.H"
#include "dummyTransform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(localPointRegion, 0);

// Reduction class to get minimum value over face.
class minEqOpFace
{
public:

    void operator()(face& x, const face& y) const
    {
        if (x.size())
        {
            label j = 0;
            forAll(x, i)
            {
                x[i] = min(x[i], y[j]);

                j = y.rcIndex(j);
            }
        }
    }
};

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Are two lists identical either in forward or in reverse order.
bool Foam::localPointRegion::isDuplicate
(
    const face& f0,
    const face& f1,
    const bool forward
)
{
    label fp1 = findIndex(f1, f0[0]);

    if (fp1 == -1)
    {
        return false;
    }

    forAll(f0, fp0)
    {
        if (f0[fp0] != f1[fp1])
        {
            return false;
        }

        if (forward)
        {
            fp1 = f1.fcIndex(fp1);
        }
        else
        {
            fp1 = f1.rcIndex(fp1);
        }
    }
    return true;
}


// Count regions per point
void Foam::localPointRegion::countPointRegions
(
    const polyMesh& mesh,
    const boolList& candidatePoint,
    const Map<label>& candidateFace,
    faceList& minRegion
)
{
    // Almost all will have only one so only
    // populate Map if more than one.
    labelList minPointRegion(mesh.nPoints(), -1);
    // From global point to local (multi-region) point numbering
    meshPointMap_.resize(candidateFace.size()/100);
    // From local (multi-region) point to regions
    DynamicList<labelList> pointRegions(meshPointMap_.size());

    // From faces with any duplicated point on it to local face
    meshFaceMap_.resize(meshPointMap_.size());

    forAllConstIter(Map<label>, candidateFace, iter)
    {
        label facei = iter.key();

        if (!mesh.isInternalFace(facei))
        {
            const face& f = mesh.faces()[facei];

            if (minRegion[facei].empty())
            {
                FatalErrorInFunction
                    << "Face from candidateFace without minRegion set." << endl
                    << "Face:" << facei << " fc:" << mesh.faceCentres()[facei]
                    << " verts:" << f << abort(FatalError);
            }

            forAll(f, fp)
            {
                label pointi = f[fp];

                // Even points which were not candidates for splitting might
                // be on multiple baffles that are being split so check.

                if (candidatePoint[pointi])
                {
                    label region = minRegion[facei][fp];

                    if (minPointRegion[pointi] == -1)
                    {
                        minPointRegion[pointi] = region;
                    }
                    else if (minPointRegion[pointi] != region)
                    {
                        // Multiple regions for this point. Add.
                        Map<label>::iterator iter = meshPointMap_.find(pointi);
                        if (iter != meshPointMap_.end())
                        {
                            labelList& regions = pointRegions[iter()];
                            if (findIndex(regions, region) == -1)
                            {
                                label sz = regions.size();
                                regions.setSize(sz+1);
                                regions[sz] = region;
                            }
                        }
                        else
                        {
                            label localPointi = meshPointMap_.size();
                            meshPointMap_.insert(pointi, localPointi);
                            labelList regions(2);
                            regions[0] = minPointRegion[pointi];
                            regions[1] = region;
                            pointRegions.append(regions);
                        }

                        label meshFaceMapI = meshFaceMap_.size();
                        meshFaceMap_.insert(facei, meshFaceMapI);
                    }
                }
            }
        }
    }
    minPointRegion.clear();

    // Add internal faces that use any duplicated point. Can only have one
    // region!
    forAllConstIter(Map<label>, candidateFace, iter)
    {
        label facei = iter.key();

        if (mesh.isInternalFace(facei))
        {
            const face& f = mesh.faces()[facei];

            forAll(f, fp)
            {
                // Note: candidatePoint test not really necessary but
                // speeds up rejection.
                if (candidatePoint[f[fp]] && meshPointMap_.found(f[fp]))
                {
                    label meshFaceMapI = meshFaceMap_.size();
                    meshFaceMap_.insert(facei, meshFaceMapI);
                }
            }
        }
    }


    // Transfer to member data
    pointRegions.shrink();
    pointRegions_.setSize(pointRegions.size());
    forAll(pointRegions, i)
    {
        pointRegions_[i].transfer(pointRegions[i]);
    }

    // Compact minRegion
    faceRegions_.setSize(meshFaceMap_.size());
    forAllConstIter(Map<label>, meshFaceMap_, iter)
    {
        faceRegions_[iter()].labelList::transfer(minRegion[iter.key()]);

        //// Print a bit
        //{
        //    label facei = iter.key();
        //    const face& f = mesh.faces()[facei];
        //    Pout<< "Face:" << facei << " fc:" << mesh.faceCentres()[facei]
        //        << " verts:" << f << endl;
        //    forAll(f, fp)
        //    {
        //        Pout<< "    " << f[fp] << " min:" << faceRegions_[iter()][fp]
        //            << endl;
        //    }
        //    Pout<< endl;
        //}
    }

    // Compact region numbering
    // ? TBD.
}


void Foam::localPointRegion::calcPointRegions
(
    const polyMesh& mesh,
    boolList& candidatePoint
)
{
    label nBnd = mesh.nFaces()-mesh.nInternalFaces();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();


    syncTools::syncPointList
    (
        mesh,
        candidatePoint,
        orEqOp<bool>(),
        false               // nullValue
    );


    // Mark any face/boundaryFace/cell with a point on a candidate point.
    // - candidateFace does not necessary have to be a baffle!
    // - candidateFace is synchronised (since candidatePoint is)
    Map<label> candidateFace(2*nBnd);
    label candidateFacei = 0;

    Map<label> candidateCell(nBnd);
    label candidateCelli = 0;

    forAll(mesh.faces(), facei)
    {
        const face& f = mesh.faces()[facei];

        forAll(f, fp)
        {
            if (candidatePoint[f[fp]])
            {
                // Mark face
                if (candidateFace.insert(facei, candidateFacei))
                {
                    candidateFacei++;
                }

                // Mark cells
                if (candidateCell.insert(faceOwner[facei], candidateCelli))
                {
                    candidateCelli++;
                }

                if (mesh.isInternalFace(facei))
                {
                    label nei = faceNeighbour[facei];
                    if (candidateCell.insert(nei, candidateCelli))
                    {
                        candidateCelli++;
                    }
                }

                break;
            }
        }
    }



    // Get global indices for cells
    globalIndex globalCells(mesh.nCells());


    // Determine for every candidate face per point the minimum region
    // (global cell) it is connected to. (candidateFaces are the
    // only ones using a
    // candidate point so the only ones that can be affected)
    faceList minRegion(mesh.nFaces());
    forAllConstIter(Map<label>, candidateFace, iter)
    {
        label facei = iter.key();
        const face& f = mesh.faces()[facei];

        if (mesh.isInternalFace(facei))
        {
            label globOwn = globalCells.toGlobal(faceOwner[facei]);
            label globNei = globalCells.toGlobal(faceNeighbour[facei]);
            minRegion[facei].setSize(f.size(), min(globOwn, globNei));
        }
        else
        {
            label globOwn = globalCells.toGlobal(faceOwner[facei]);
            minRegion[facei].setSize(f.size(), globOwn);
        }
    }

    // Now minimise over all faces that are connected through internal
    // faces or through cells. This loop iterates over the max number of
    // cells connected to a point (=8 for hex mesh)

    while (true)
    {
        // Transport minimum from face across cell
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Map<label> minPointValue(100);
        label nChanged = 0;
        forAllConstIter(Map<label>, candidateCell, iter)
        {
            minPointValue.clear();

            label celli = iter.key();
            const cell& cFaces = mesh.cells()[celli];

            // Determine minimum per point
            forAll(cFaces, cFacei)
            {
                label facei = cFaces[cFacei];

                if (minRegion[facei].size())
                {
                    const face& f = mesh.faces()[facei];

                    forAll(f, fp)
                    {
                        label pointi = f[fp];
                        Map<label>::iterator iter = minPointValue.find(pointi);

                        if (iter == minPointValue.end())
                        {
                            minPointValue.insert(pointi, minRegion[facei][fp]);
                        }
                        else
                        {
                            label currentMin = iter();
                            iter() = min(currentMin, minRegion[facei][fp]);
                        }
                    }
                }
            }

            // Set face minimum from point minimum
            forAll(cFaces, cFacei)
            {
                label facei = cFaces[cFacei];

                if (minRegion[facei].size())
                {
                    const face& f = mesh.faces()[facei];

                    forAll(f, fp)
                    {
                        label minVal = minPointValue[f[fp]];

                        if (minVal != minRegion[facei][fp])
                        {
                            minRegion[facei][fp] = minVal;
                            nChanged++;
                        }
                    }
                }
            }
        }

        // Pout<< "nChanged:" << nChanged << endl;

        if (returnReduce(nChanged, sumOp<label>()) == 0)
        {
            break;
        }


        // Transport minimum across coupled faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        SubList<face> l
        (
            minRegion,
            mesh.nFaces()-mesh.nInternalFaces(),
            mesh.nInternalFaces()
        );
        syncTools::syncBoundaryFaceList
        (
            mesh,
            l,
            minEqOpFace(),
            Foam::dummyTransform()  // dummy transformation
        );
    }


    // Count regions per point
    countPointRegions(mesh, candidatePoint, candidateFace, minRegion);
    minRegion.clear();


    //// Print points with multiple regions. These points need to be duplicated.
    // forAllConstIter(Map<label>, meshPointMap_, iter)
    //{
    //    Pout<< "point:" << iter.key()
    //        << " coord:" << mesh.points()[iter.key()]
    //        << " regions:" << pointRegions_[iter()] << endl;
    //}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::localPointRegion::localPointRegion(const polyMesh& mesh)
:
    // nRegions_(0),
    meshPointMap_(0),
    pointRegions_(0),
    meshFaceMap_(0),
    faceRegions_(0)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Get any point on the outside which is on a non-coupled boundary
    boolList candidatePoint(mesh.nPoints(), false);

    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {
            const polyPatch& pp = patches[patchi];

            forAll(pp.meshPoints(), i)
            {
                candidatePoint[pp.meshPoints()[i]] = true;
            }
        }
    }

    calcPointRegions(mesh, candidatePoint);
}


Foam::localPointRegion::localPointRegion
(
    const polyMesh& mesh,
    const labelList& candidatePoints
)
:
    // nRegions_(0),
    meshPointMap_(0),
    pointRegions_(0),
    meshFaceMap_(0),
    faceRegions_(0)
{
    boolList candidatePoint(mesh.nPoints(), false);

    forAll(candidatePoints, i)
    {
        candidatePoint[candidatePoints[i]] = true;
    }

    calcPointRegions(mesh, candidatePoint);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return a list (in allPatch indices) with either -1 or the face label
// of the face that uses the same vertices.
Foam::labelList Foam::localPointRegion::findDuplicateFaces
(
    const primitiveMesh& mesh,
    const labelList& boundaryFaces
)
{
    // Addressing engine for all boundary faces.
    indirectPrimitivePatch allPatch
    (
        IndirectList<face>(mesh.faces(), boundaryFaces),
        mesh.points()
    );

    labelList duplicateFace(allPatch.size(), -1);
    label nDuplicateFaces = 0;

    // Find all duplicate faces.
    forAll(allPatch, bFacei)
    {
        const face& f = allPatch.localFaces()[bFacei];

        // Get faces connected to f[0].
        // Check whether share all points with f.
        const labelList& pFaces = allPatch.pointFaces()[f[0]];

        forAll(pFaces, i)
        {
            label otherFacei = pFaces[i];

            if (otherFacei > bFacei)
            {
                const face& otherF = allPatch.localFaces()[otherFacei];

                if (isDuplicate(f, otherF, true))
                {
                    FatalErrorInFunction
                        << "Face:" << bFacei + mesh.nInternalFaces()
                        << " has local points:" << f
                        << " which are in same order as face:"
                        << otherFacei + mesh.nInternalFaces()
                        << " with local points:" << otherF
                        << abort(FatalError);
                }
                else if (isDuplicate(f, otherF, false))
                {
                    label meshFace0 = bFacei + mesh.nInternalFaces();
                    label meshFace1 = otherFacei + mesh.nInternalFaces();

                    if
                    (
                        duplicateFace[bFacei] != -1
                     || duplicateFace[otherFacei] != -1
                    )
                    {
                        FatalErrorInFunction
                            << "One of two duplicate faces already marked"
                            << " as duplicate." << nl
                            << "This means that three or more faces share"
                            << " the same points and this is illegal." << nl
                            << "Face:" << meshFace0
                            << " with local points:" << f
                            << " which are in same order as face:"
                            << meshFace1
                            << " with local points:" << otherF
                            << abort(FatalError);
                    }

                    duplicateFace[bFacei] = otherFacei;
                    duplicateFace[otherFacei] = bFacei;
                    nDuplicateFaces++;
                }
            }
        }
    }

    return duplicateFace;
}


Foam::List<Foam::labelPair> Foam::localPointRegion::findDuplicateFacePairs
(
    const polyMesh& mesh
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Faces to test: all boundary faces
    labelList testFaces
    (
        identity(mesh.nFaces()-mesh.nInternalFaces())
      + mesh.nInternalFaces()
    );

    // Find corresponding baffle face (or -1)
    const labelList duplicateFace(findDuplicateFaces(mesh, testFaces));

    // Convert into list of coupled face pairs (mesh face labels).
    DynamicList<labelPair> baffles(testFaces.size());

    forAll(duplicateFace, i)
    {
        label otherFacei = duplicateFace[i];

        if (otherFacei != -1 && i < otherFacei)
        {
            label meshFace0 = testFaces[i];
            label patch0 = patches.whichPatch(meshFace0);
            label meshFace1 = testFaces[otherFacei];
            label patch1 = patches.whichPatch(meshFace1);

            // Check for illegal topology. Should normally not happen!
            if
            (
                (patch0 != -1 && isA<processorPolyPatch>(patches[patch0]))
             || (patch1 != -1 && isA<processorPolyPatch>(patches[patch1]))
            )
            {
                FatalErrorInFunction
                    << "One of two duplicate faces is on"
                    << " processorPolyPatch."
                    << "This is not allowed." << nl
                    << "Face:" << meshFace0
                    << " is on patch:" << patches[patch0].name()
                    << nl
                    << "Face:" << meshFace1
                    << " is on patch:" << patches[patch1].name()
                    << abort(FatalError);
            }

            baffles.append(labelPair(meshFace0, meshFace1));
        }
    }
    return baffles.shrink();
}


void Foam::localPointRegion::topoChange(const polyTopoChangeMap& map)
{
    {
        Map<label> newMap(meshFaceMap_.size());

        forAllConstIter(Map<label>, meshFaceMap_, iter)
        {
            label newFacei = map.reverseFaceMap()[iter.key()];

            if (newFacei >= 0)
            {
                newMap.insert(newFacei, iter());
            }
        }
        meshFaceMap_.transfer(newMap);
    }
    {
        Map<label> newMap(meshPointMap_.size());

        forAllConstIter(Map<label>, meshPointMap_, iter)
        {
            label newPointi = map.reversePointMap()[iter.key()];

            if (newPointi >= 0)
            {
                newMap.insert(newPointi, iter());
            }
        }

        meshPointMap_.transfer(newMap);
    }
}


// ************************************************************************* //
