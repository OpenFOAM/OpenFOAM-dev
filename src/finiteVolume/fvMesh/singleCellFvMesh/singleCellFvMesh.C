/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "singleCellFvMesh.H"
#include "syncTools.H"
#include "uindirectPrimitivePatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::singleCellFvMesh::agglomerateMesh
(
    const fvMesh& mesh,
    const labelListList& agglom
)
{
    // Conversion is a two step process:
    // - from original (fine) patch faces to agglomerations (aggloms might not
    //   be in correct patch order)
    // - from agglomerations to coarse patch faces

    const polyBoundaryMesh& oldPatches = mesh.boundaryMesh();

    // Check agglomeration within patch face range and continuous
    labelList nAgglom(oldPatches.size(), 0);

    forAll(oldPatches, patchi)
    {
        const polyPatch& pp = oldPatches[patchi];
        if (pp.size() > 0)
        {
            nAgglom[patchi] = max(agglom[patchi])+1;

            forAll(pp, i)
            {
                if (agglom[patchi][i] < 0  || agglom[patchi][i] >= pp.size())
                {
                    FatalErrorInFunction
                        << "agglomeration on patch " << patchi
                        << " is out of range 0.." << pp.size()-1
                        << exit(FatalError);
                }
            }
        }
    }

    // Check agglomeration is sync
    {
        // Get neighbouring agglomeration
        labelList nbrAgglom(mesh.nFaces()-mesh.nInternalFaces());
        forAll(oldPatches, patchi)
        {
            const polyPatch& pp = oldPatches[patchi];

            if (pp.coupled())
            {
                label offset = pp.start()-mesh.nInternalFaces();
                forAll(pp, i)
                {
                    nbrAgglom[offset+i] = agglom[patchi][i];
                }
            }
        }
        syncTools::swapBoundaryFaceList(mesh, nbrAgglom);


        // Get correspondence between this agglomeration and remote one
        Map<label> localToNbr(nbrAgglom.size()/10);

        forAll(oldPatches, patchi)
        {
            const polyPatch& pp = oldPatches[patchi];

            if (pp.coupled())
            {
                label offset = pp.start()-mesh.nInternalFaces();

                forAll(pp, i)
                {
                    label bFacei = offset+i;
                    label myZone = agglom[patchi][i];
                    label nbrZone = nbrAgglom[bFacei];

                    Map<label>::const_iterator iter = localToNbr.find(myZone);

                    if (iter == localToNbr.end())
                    {
                        // First occurrence of this zone. Store correspondence
                        // to remote zone number.
                        localToNbr.insert(myZone, nbrZone);
                    }
                    else
                    {
                        // Check that zone numbers are still the same.
                        if (iter() != nbrZone)
                        {
                            FatalErrorInFunction
                                << "agglomeration is not synchronised across"
                                << " coupled patch " << pp.name()
                                << endl
                                << "Local agglomeration " << myZone
                                << ". Remote agglomeration " << nbrZone
                                << exit(FatalError);
                        }
                    }
                }
            }
        }
    }


    label coarseI = 0;
    forAll(nAgglom, patchi)
    {
        coarseI += nAgglom[patchi];
    }
    // New faces
    faceList patchFaces(coarseI);
    // New patch start and size
    labelList patchStarts(oldPatches.size());
    labelList patchSizes(oldPatches.size());

    // From new patch face back to agglomeration
    patchFaceMap_.setSize(oldPatches.size());

    // From fine face to coarse face (or -1)
    reverseFaceMap_.setSize(mesh.nFaces());
    reverseFaceMap_.labelList::operator=(-1);

    // Face counter
    coarseI = 0;


    forAll(oldPatches, patchi)
    {
        patchStarts[patchi] = coarseI;

        const polyPatch& pp = oldPatches[patchi];

        if (pp.size() > 0)
        {
            patchFaceMap_[patchi].setSize(nAgglom[patchi]);

            // Patchfaces per agglomeration
            labelListList agglomToPatch
            (
                invertOneToMany(nAgglom[patchi], agglom[patchi])
            );

            // From agglomeration to compact patch face
            labelList agglomToFace(nAgglom[patchi], -1);

            forAll(pp, i)
            {
                label myAgglom = agglom[patchi][i];

                if (agglomToFace[myAgglom] == -1)
                {
                    // Agglomeration not yet done. We now have:
                    // - coarseI                  : current coarse mesh face
                    // - patchStarts[patchi]      : coarse mesh patch start
                    // - myAgglom                 : agglomeration
                    // -  agglomToPatch[myAgglom] : fine mesh faces for zone
                    label coarsePatchFacei = coarseI - patchStarts[patchi];
                    patchFaceMap_[patchi][coarsePatchFacei] = myAgglom;
                    agglomToFace[myAgglom] = coarsePatchFacei;

                    const labelList& fineFaces = agglomToPatch[myAgglom];

                    // Create overall map from fine mesh faces to coarseI.
                    forAll(fineFaces, fineI)
                    {
                        reverseFaceMap_[pp.start()+fineFaces[fineI]] = coarseI;
                    }

                    // Construct single face
                    uindirectPrimitivePatch upp
                    (
                        UIndirectList<face>(pp, fineFaces),
                        pp.points()
                    );

                    if (upp.edgeLoops().size() != 1)
                    {
                        FatalErrorInFunction
                            << "agglomeration does not create a"
                            << " single, non-manifold"
                            << " face for agglomeration " << myAgglom
                            << " on patch " <<  patchi
                            << exit(FatalError);
                    }

                    patchFaces[coarseI++] = face
                    (
                        renumber
                        (
                            upp.meshPoints(),
                            upp.edgeLoops()[0]
                        )
                    );
                }
            }
        }

        patchSizes[patchi] = coarseI-patchStarts[patchi];
    }

    // Pout<< "patchStarts:" << patchStarts << endl;
    // Pout<< "patchSizes:" << patchSizes << endl;

    // Compact numbering for points
    reversePointMap_.setSize(mesh.nPoints());
    reversePointMap_.labelList::operator=(-1);
    label newI = 0;

    forAll(patchFaces, coarseI)
    {
        face& f = patchFaces[coarseI];

        forAll(f, fp)
        {
            if (reversePointMap_[f[fp]] == -1)
            {
                reversePointMap_[f[fp]] = newI++;
            }

            f[fp] = reversePointMap_[f[fp]];
        }
    }

    pointMap_ = invert(newI, reversePointMap_);

    // Subset used points
    pointField boundaryPoints(mesh.points(), pointMap_);

    // Add patches (on still zero sized mesh)
    List<polyPatch*> newPatches(oldPatches.size());
    forAll(oldPatches, patchi)
    {
        newPatches[patchi] = oldPatches[patchi].clone
        (
            boundaryMesh(),
            patchi,
            0,
            0
        ).ptr();
    }
    addFvPatches(newPatches);

    // Owner, neighbour is trivial
    labelList owner(patchFaces.size(), 0);
    labelList neighbour(0);


    // actually change the mesh
    resetPrimitives
    (
        std::move(boundaryPoints),
        std::move(patchFaces),
        std::move(owner),
        std::move(neighbour),
        patchSizes,
        patchStarts,
        true                // syncPar
    );


    // Adapt the zones
    cellZones().clear();
    cellZones().setSize(mesh.cellZones().size());
    {
        forAll(mesh.cellZones(), zoneI)
        {
            const cellZone& oldFz = mesh.cellZones()[zoneI];

            DynamicList<label> newAddressing;

            // Note: uncomment if you think it makes sense. Note that value
            // of cell0 is the average.
            //// Was old cell0 in this cellZone?
            // if (oldFz.localID(0) != -1)
            //{
            //    newAddressing.append(0);
            //}

            cellZones().set
            (
                zoneI,
                oldFz.clone
                (
                    newAddressing,
                    zoneI,
                    cellZones()
                )
            );
        }
    }

    faceZones().clear();
    faceZones().setSize(mesh.faceZones().size());
    {
        forAll(mesh.faceZones(), zoneI)
        {
            const faceZone& oldFz = mesh.faceZones()[zoneI];

            DynamicList<label> newAddressing(oldFz.size());
            DynamicList<bool> newFlipMap(oldFz.size());

            forAll(oldFz, i)
            {
                label newFacei = reverseFaceMap_[oldFz[i]];

                if (newFacei != -1)
                {
                    newAddressing.append(newFacei);
                    newFlipMap.append(oldFz.flipMap()[i]);
                }
            }

            faceZones().set
            (
                zoneI,
                oldFz.clone
                (
                    newAddressing,
                    newFlipMap,
                    zoneI,
                    faceZones()
                )
            );
        }
    }


    pointZones().clear();
    pointZones().setSize(mesh.pointZones().size());
    {
        forAll(mesh.pointZones(), zoneI)
        {
            const pointZone& oldFz = mesh.pointZones()[zoneI];

            DynamicList<label> newAddressing(oldFz.size());

            forAll(oldFz, i)
            {
                label newPointi  = reversePointMap_[oldFz[i]];
                if (newPointi != -1)
                {
                    newAddressing.append(newPointi);
                }
            }

            pointZones().set
            (
                zoneI,
                oldFz.clone
                (
                    pointZones(),
                    zoneI,
                    newAddressing
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singleCellFvMesh::singleCellFvMesh
(
    const IOobject& io,
    const fvMesh& mesh
)
:
    fvMesh
    (
        io,
        pointField(), // points
        faceList(),   // faces
        labelList(),  // allOwner
        labelList(),  // allNeighbour
        false         // syncPar
    ),
    patchFaceAgglomeration_
    (
        IOobject
        (
            "patchFaceAgglomeration",
            io.instance(),
            fvMesh::meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        ),
        0
    ),
    patchFaceMap_
    (
        IOobject
        (
            "patchFaceMap",
            io.instance(),
            fvMesh::meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        ),
        mesh.boundaryMesh().size()
    ),
    reverseFaceMap_
    (
        IOobject
        (
            "reverseFaceMap",
            io.instance(),
            fvMesh::meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        ),
        mesh.nFaces()
    ),
    pointMap_
    (
        IOobject
        (
            "pointMap",
            io.instance(),
            fvMesh::meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        ),
        mesh.nPoints()
    ),
    reversePointMap_
    (
        IOobject
        (
            "reversePointMap",
            io.instance(),
            fvMesh::meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        ),
        mesh.nPoints()
    )
{
    const polyBoundaryMesh& oldPatches = mesh.boundaryMesh();

    labelListList agglom(oldPatches.size());

    forAll(oldPatches, patchi)
    {
        agglom[patchi] = identityMap(oldPatches[patchi].size());
    }

    agglomerateMesh(mesh, agglom);
}


Foam::singleCellFvMesh::singleCellFvMesh
(
    const IOobject& io,
    const fvMesh& mesh,
    const labelListList& patchFaceAgglomeration
)
:
    fvMesh
    (
        io,
        pointField(), // points
        faceList(),   // faces
        labelList(),  // allOwner
        labelList(),  // allNeighbour
        false         // syncPar
    ),
    patchFaceAgglomeration_
    (
        IOobject
        (
            "patchFaceAgglomeration",
            io.instance(),
            fvMesh::meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        ),
        patchFaceAgglomeration
    ),
    patchFaceMap_
    (
        IOobject
        (
            "patchFaceMap",
            io.instance(),
            fvMesh::meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        ),
        mesh.boundaryMesh().size()
    ),
    reverseFaceMap_
    (
        IOobject
        (
            "reverseFaceMap",
            io.instance(),
            fvMesh::meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        ),
        mesh.nFaces()
    ),
    pointMap_
    (
        IOobject
        (
            "pointMap",
            io.instance(),
            fvMesh::meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        ),
        mesh.nPoints()
    ),
    reversePointMap_
    (
        IOobject
        (
            "reversePointMap",
            io.instance(),
            fvMesh::meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        ),
        mesh.nPoints()
    )
{
    agglomerateMesh(mesh, patchFaceAgglomeration);
}


Foam::singleCellFvMesh::singleCellFvMesh(const IOobject& io)
:
    fvMesh(io),
    patchFaceAgglomeration_
    (
        IOobject
        (
            "patchFaceAgglomeration",
            io.instance(),
            fvMesh::meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        )
    ),
    patchFaceMap_
    (
        IOobject
        (
            "patchFaceMap",
            io.instance(),
            fvMesh::meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        )
    ),
    reverseFaceMap_
    (
        IOobject
        (
            "reverseFaceMap",
            io.instance(),
            fvMesh::meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        )
    ),
    pointMap_
    (
        IOobject
        (
            "pointMap",
            io.instance(),
            fvMesh::meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        )
    ),
    reversePointMap_
    (
        IOobject
        (
            "reversePointMap",
            io.instance(),
            fvMesh::meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        )
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// ************************************************************************* //
