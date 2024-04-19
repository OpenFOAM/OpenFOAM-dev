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

\*---------------------------------------------------------------------------*/

#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "processorPolyPatch.H"
#include "stringListOps.H"
#include "PstreamBuffers.H"
#include "lduSchedule.H"
#include "globalMeshData.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyBoundaryMesh, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyBoundaryMesh::polyBoundaryMesh
(
    const IOobject& io,
    const polyMesh& mesh
)
:
    polyPatchList(),
    regIOobject(io),
    mesh_(mesh)
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
    )
    {
        if (readOpt() == IOobject::MUST_READ_IF_MODIFIED)
        {
            WarningInFunction
                << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
                << " does not support automatic rereading."
                << endl;
        }


        polyPatchList& patches = *this;

        // Read polyPatchList
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        patches.setSize(patchEntries.size());

        forAll(patches, patchi)
        {
            patches.set
            (
                patchi,
                polyPatch::New
                (
                    patchEntries[patchi].keyword(),
                    patchEntries[patchi].dict(),
                    patchi,
                    *this
                )
            );
        }

        // Check state of IOstream
        is.check
        (
            "polyBoundaryMesh::polyBoundaryMesh"
            "(const IOobject&, const polyMesh&)"
        );

        close();
    }
}


Foam::polyBoundaryMesh::polyBoundaryMesh
(
    const IOobject& io,
    const polyMesh& pm,
    const label size
)
:
    polyPatchList(size),
    regIOobject(io),
    mesh_(pm)
{}


Foam::polyBoundaryMesh::polyBoundaryMesh
(
    const IOobject& io,
    const polyMesh& pm,
    const polyPatchList& ppl
)
:
    polyPatchList(),
    regIOobject(io),
    mesh_(pm)
{
    if
    (
        (this->readOpt() == IOobject::READ_IF_PRESENT && this->headerOk())
     || this->readOpt() == IOobject::MUST_READ
     || this->readOpt() == IOobject::MUST_READ_IF_MODIFIED
    )
    {

        if (readOpt() == IOobject::MUST_READ_IF_MODIFIED)
        {
            WarningInFunction
                << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
                << " does not support automatic rereading."
                << endl;
        }

        polyPatchList& patches = *this;

        // Read polyPatchList
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        patches.setSize(patchEntries.size());

        forAll(patches, patchi)
        {
            patches.set
            (
                patchi,
                polyPatch::New
                (
                    patchEntries[patchi].keyword(),
                    patchEntries[patchi].dict(),
                    patchi,
                    *this
                )
            );
        }

        // Check state of IOstream
        is.check
        (
            "polyBoundaryMesh::polyBoundaryMesh"
            "(const IOobject&, const polyMesh&, const polyPatchList&)"
        );

        close();
    }
    else
    {
        polyPatchList& patches = *this;
        patches.setSize(ppl.size());
        forAll(patches, patchi)
        {
            patches.set(patchi, ppl[patchi].clone(*this).ptr());
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyBoundaryMesh::~polyBoundaryMesh()
{}


void Foam::polyBoundaryMesh::clearGeom()
{
    forAll(*this, patchi)
    {
        if (this->set(patchi))
        {
            operator[](patchi).clearGeom();
        }
    }
}


void Foam::polyBoundaryMesh::clearAddressing()
{
    nbrEdgesPtr_.clear();
    patchIndicesPtr_.clear();
    patchFaceIndicesPtr_.clear();
    groupPatchIndicesPtr_.clear();

    forAll(*this, patchi)
    {
        if (this->set(patchi))
        {
            operator[](patchi).clearAddressing();
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyBoundaryMesh::calcGeometry()
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchi)
        {
            operator[](patchi).initCalcGeometry(pBufs);
        }

        pBufs.finishedSends();

        forAll(*this, patchi)
        {
            operator[](patchi).calcGeometry(pBufs);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        forAll(patchSchedule, patchEvali)
        {
            const label patchi = patchSchedule[patchEvali].patch;

            if (patchSchedule[patchEvali].init)
            {
                operator[](patchi).initCalcGeometry(pBufs);
            }
            else
            {
                operator[](patchi).calcGeometry(pBufs);
            }
        }
    }
}


const Foam::List<Foam::labelPairList>&
Foam::polyBoundaryMesh::nbrEdges() const
{
    if (Pstream::parRun())
    {
        WarningInFunction
            << "Neighbour edge addressing not correct across parallel"
            << " boundaries." << endl;
    }

    if (!nbrEdgesPtr_.valid())
    {
        nbrEdgesPtr_.reset(new List<labelPairList>(size()));
        List<labelPairList>& nbrEdges = nbrEdgesPtr_();

        // Initialise.
        label nEdgePairs = 0;
        forAll(*this, patchi)
        {
            const polyPatch& pp = operator[](patchi);

            nbrEdges[patchi].setSize(pp.nEdges() - pp.nInternalEdges());

            forAll(nbrEdges[patchi], i)
            {
                labelPair& edgeInfo = nbrEdges[patchi][i];

                edgeInfo[0] = -1;
                edgeInfo[1] = -1;
            }

            nEdgePairs += pp.nEdges() - pp.nInternalEdges();
        }

        // From mesh edge (expressed as a point pair so as not to construct
        // point addressing) to patch + relative edge index.
        HashTable<labelPair, edge, Hash<edge>> pointsToEdge(nEdgePairs);

        forAll(*this, patchi)
        {
            const polyPatch& pp = operator[](patchi);

            const edgeList& edges = pp.edges();

            for
            (
                label edgei = pp.nInternalEdges();
                edgei < edges.size();
                edgei++
            )
            {
                // Edge in patch local points
                const edge& e = edges[edgei];

                // Edge in mesh points.
                edge meshEdge(pp.meshPoints()[e[0]], pp.meshPoints()[e[1]]);

                HashTable<labelPair, edge, Hash<edge>>::iterator fnd =
                    pointsToEdge.find(meshEdge);

                if (fnd == pointsToEdge.end())
                {
                    // First occurrence of mesh edge. Store patch and my
                    // local index.
                    pointsToEdge.insert
                    (
                        meshEdge,
                        labelPair
                        (
                            patchi,
                            edgei - pp.nInternalEdges()
                        )
                    );
                }
                else
                {
                    // Second occurrence. Store.
                    const labelPair& edgeInfo = fnd();

                    nbrEdges[patchi][edgei - pp.nInternalEdges()] =
                        edgeInfo;

                    nbrEdges[edgeInfo[0]][edgeInfo[1]]
                         = labelPair(patchi, edgei - pp.nInternalEdges());

                    // Found all two occurrences of this edge so remove from
                    // hash to save space. Note that this will give lots of
                    // problems if the polyBoundaryMesh is multiply connected.
                    pointsToEdge.erase(meshEdge);
                }
            }
        }

        if (pointsToEdge.size())
        {
            FatalErrorInFunction
                << "Not all boundary edges of patches match up." << nl
                << "Is the outside of your mesh multiply connected?"
                << abort(FatalError);
        }

        forAll(*this, patchi)
        {
            const polyPatch& pp = operator[](patchi);

            const labelPairList& nbrEdgesp = nbrEdges[patchi];

            forAll(nbrEdgesp, i)
            {
                const labelPair& edgeInfo = nbrEdgesp[i];

                if (edgeInfo[0] == -1 || edgeInfo[1] == -1)
                {
                    label edgeI = pp.nInternalEdges() + i;
                    const edge& e = pp.edges()[edgeI];

                    FatalErrorInFunction
                        << "Not all boundary edges of patches match up." << nl
                        << "Edge " << edgeI << " on patch " << pp.name()
                        << " end points " << pp.localPoints()[e[0]] << ' '
                        << pp.localPoints()[e[1]] << " is not matched to an"
                        << " edge on any other patch." << nl
                        << "Is the outside of your mesh multiply connected?"
                        << abort(FatalError);
                }
            }
        }
    }

    return nbrEdgesPtr_();
}


const Foam::labelList& Foam::polyBoundaryMesh::patchIndices() const
{
    if (!patchIndicesPtr_.valid())
    {
        patchIndicesPtr_.reset
        (
            new labelList
            (
                mesh_.nFaces()
              - mesh_.nInternalFaces()
            )
        );
        labelList& patchIndices = patchIndicesPtr_();

        const polyBoundaryMesh& bm = *this;

        forAll(bm, patchi)
        {
            label bFacei = bm[patchi].start() - mesh_.nInternalFaces();
            forAll(bm[patchi], i)
            {
                patchIndices[bFacei++] = patchi;
            }
        }
    }
    return patchIndicesPtr_();
}


const Foam::labelList& Foam::polyBoundaryMesh::patchFaceIndices() const
{
    if (!patchFaceIndicesPtr_.valid())
    {
        patchFaceIndicesPtr_.reset
        (
            new labelList
            (
                mesh_.nFaces()
              - mesh_.nInternalFaces()
            )
        );
        labelList& patchFaceID = patchFaceIndicesPtr_();

        const polyBoundaryMesh& bm = *this;

        forAll(bm, patchi)
        {
            label bFacei = bm[patchi].start() - mesh_.nInternalFaces();
            forAll(bm[patchi], i)
            {
                patchFaceID[bFacei++] = i;
            }
        }
    }
    return patchFaceIndicesPtr_();
}


const Foam::HashTable<Foam::labelList, Foam::word>&
Foam::polyBoundaryMesh::groupPatchIndices() const
{
    if (!groupPatchIndicesPtr_.valid())
    {
        groupPatchIndicesPtr_.reset(new HashTable<labelList, word>(10));
        HashTable<labelList, word>& groupPatchIndices = groupPatchIndicesPtr_();

        const polyBoundaryMesh& bm = *this;

        forAll(bm, patchi)
        {
            const wordList& groups = bm[patchi].inGroups();

            forAll(groups, i)
            {
                const word& name = groups[i];

                HashTable<labelList, word>::iterator iter
                (
                    groupPatchIndices.find
                    (
                        name
                    )
                );

                if (iter != groupPatchIndices.end())
                {
                    iter().append(patchi);
                }
                else
                {
                    groupPatchIndices.insert(name, labelList(1, patchi));
                }
            }
        }

        // Remove patch names from patchGroups
        forAll(bm, patchi)
        {
            if (groupPatchIndices.erase(bm[patchi].name()))
            {
                WarningInFunction
                    << "Removing patchGroup '" << bm[patchi].name()
                    << "' which clashes with patch " << patchi
                    << " of the same name."
                    << endl;
            }
        }
    }

    return groupPatchIndicesPtr_();
}


void Foam::polyBoundaryMesh::setGroup
(
    const word& groupName,
    const labelList& patchIndices
)
{
    groupPatchIndicesPtr_.clear();

    polyPatchList& patches = *this;

    boolList donePatch(patches.size(), false);

    // Add to specified patches
    forAll(patchIndices, i)
    {
        label patchi = patchIndices[i];
        polyPatch& pp = patches[patchi];

        if (!pp.inGroup(groupName))
        {
            pp.inGroups().append(groupName);
        }
        donePatch[patchi] = true;
    }

    // Remove from other patches
    forAll(patches, patchi)
    {
        if (!donePatch[patchi])
        {
            polyPatch& pp = patches[patchi];

            label newI = 0;
            if (pp.inGroup(groupName))
            {
                wordList& groups = pp.inGroups();

                forAll(groups, i)
                {
                    if (groups[i] != groupName)
                    {
                        groups[newI++] = groups[i];
                    }
                }
                groups.setSize(newI);
            }
        }
    }
}


Foam::wordList Foam::polyBoundaryMesh::toc() const
{
    const polyPatchList& patches = *this;

    wordList t(patches.size());

    forAll(patches, patchi)
    {
        t[patchi] = patches[patchi].name();
    }

    return t;
}


Foam::wordList Foam::polyBoundaryMesh::names() const
{
    return toc();
}


Foam::wordList Foam::polyBoundaryMesh::types() const
{
    const polyPatchList& patches = *this;

    wordList t(patches.size());

    forAll(patches, patchi)
    {
        t[patchi] = patches[patchi].type();
    }

    return t;
}


Foam::wordList Foam::polyBoundaryMesh::physicalTypes() const
{
    const polyPatchList& patches = *this;

    wordList t(patches.size());

    forAll(patches, patchi)
    {
        t[patchi] = patches[patchi].physicalType();
    }

    return t;
}


Foam::labelList Foam::polyBoundaryMesh::findIndices
(
    const wordRe& key,
    const bool usePatchGroups
) const
{
    DynamicList<label> indices;

    if (!key.empty())
    {
        if (key.isPattern())
        {
            indices = findStrings(key, this->names());

            if (usePatchGroups && groupPatchIndices().size())
            {
                labelHashSet indexSet(indices);

                const wordList allGroupNames = groupPatchIndices().toc();
                labelList groupIndices = findStrings(key, allGroupNames);
                forAll(groupIndices, i)
                {
                    const word& grpName = allGroupNames[groupIndices[i]];
                    const labelList& patchIndices =
                        groupPatchIndices()[grpName];
                    forAll(patchIndices, j)
                    {
                        if (indexSet.insert(patchIndices[j]))
                        {
                            indices.append(patchIndices[j]);
                        }
                    }
                }
            }
        }
        else
        {
            // Literal string. Special version of above to avoid
            // unnecessary memory allocations

            indices.setCapacity(1);
            forAll(*this, i)
            {
                if (key == operator[](i).name())
                {
                    indices.append(i);
                    break;
                }
            }

            if (usePatchGroups && groupPatchIndices().size())
            {
                const HashTable<labelList, word>::const_iterator iter =
                    groupPatchIndices().find(key);

                if (iter != groupPatchIndices().end())
                {
                    labelHashSet indexSet(indices);

                    const labelList& patchIndices = iter();
                    forAll(patchIndices, j)
                    {
                        if (indexSet.insert(patchIndices[j]))
                        {
                            indices.append(patchIndices[j]);
                        }
                    }
                }
            }
        }
    }

    return indices;
}


Foam::label Foam::polyBoundaryMesh::findIndex(const word& patchName) const
{
    const polyPatchList& patches = *this;

    forAll(patches, patchi)
    {
        if (patches[patchi].name() == patchName)
        {
            return patchi;
        }
    }

    // Patch not found
    if (debug)
    {
        Pout<< "label polyBoundaryMesh::findIndex(const word&) const"
            << "Patch named " << patchName << " not found.  "
            << "List of available patch names: " << names() << endl;
    }

    // Not found, return -1
    return -1;
}


Foam::label Foam::polyBoundaryMesh::whichPatch(const label faceIndex) const
{
    // Find out which patch the current face belongs to by comparing label
    // with patch start labels.
    // If the face is internal, return -1;
    // if it is off the end of the list, abort
    if (faceIndex < mesh().nInternalFaces())
    {
        return -1;
    }
    else if (faceIndex >= mesh().nFaces())
    {
        FatalErrorInFunction
            << "given label " << faceIndex
            << " greater than the number of geometric faces " << mesh().nFaces()
            << abort(FatalError);
    }


    forAll(*this, patchi)
    {
        const polyPatch& bp = operator[](patchi);

        if
        (
            faceIndex >= bp.start()
         && faceIndex < bp.start() + bp.size()
        )
        {
            return patchi;
        }
    }

    // If not in any of above, it is trouble!
    FatalErrorInFunction
        << "Cannot find face " << faceIndex << " in any of the patches "
        << names() << nl
        << "It seems your patches are not consistent with the mesh :"
        << " internalFaces:" << mesh().nInternalFaces()
        << "  total number of faces:" << mesh().nFaces()
        << abort(FatalError);

    return -1;
}


Foam::labelHashSet Foam::polyBoundaryMesh::patchSet
(
    const UList<wordRe>& patchNames,
    const bool warnNotFound,
    const bool usePatchGroups
) const
{
    const wordList allPatchNames(this->names());
    labelHashSet ids(size());

    forAll(patchNames, i)
    {
        const wordRe& patchName = patchNames[i];

        // Treat the given patch names as wild-cards and search the set
        // of all patch names for matches
        labelList patchIndices = findStrings(patchName, allPatchNames);

        forAll(patchIndices, j)
        {
            ids.insert(patchIndices[j]);
        }

        if (patchIndices.empty())
        {
            if (usePatchGroups)
            {
                const wordList allGroupNames = groupPatchIndices().toc();

                // Regard as group name
                labelList groupIndices = findStrings(patchName, allGroupNames);

                forAll(groupIndices, i)
                {
                    const word& name = allGroupNames[groupIndices[i]];
                    const labelList& extraPatchIndices =
                        groupPatchIndices()[name];

                    forAll(extraPatchIndices, extraI)
                    {
                        ids.insert(extraPatchIndices[extraI]);
                    }
                }

                if (groupIndices.empty() && warnNotFound)
                {
                    WarningInFunction
                        << "Cannot find any patch or group names matching "
                        << patchName
                        << endl;
                }
            }
            else if (warnNotFound)
            {
                WarningInFunction
                    << "Cannot find any patch names matching " << patchName
                    << endl;
            }
        }
    }

    return ids;
}


void Foam::polyBoundaryMesh::matchGroups
(
    const labelUList& patchIndices,
    wordList& groups,
    labelHashSet& nonGroupPatches
) const
{
    // Current matched groups
    DynamicList<word> matchedGroups(1);

    // Current set of unmatched patches
    nonGroupPatches = labelHashSet(patchIndices);

    const HashTable<labelList, word>& groupPatchIndices =
        this->groupPatchIndices();

    for
    (
        HashTable<labelList,word>::const_iterator iter =
            groupPatchIndices.begin();
        iter != groupPatchIndices.end();
        ++iter
    )
    {
        // Store currently unmatched patches so we can restore
        labelHashSet oldNonGroupPatches(nonGroupPatches);

        // Match by deleting patches in group from the current set and seeing
        // if all have been deleted.
        labelHashSet groupPatchSet(iter());

        label nMatch = nonGroupPatches.erase(groupPatchSet);

        if (nMatch == groupPatchSet.size())
        {
            matchedGroups.append(iter.key());
        }
        else if (nMatch != 0)
        {
            // No full match. Undo.
            nonGroupPatches.transfer(oldNonGroupPatches);
        }
    }

    groups.transfer(matchedGroups);
}


bool Foam::polyBoundaryMesh::checkParallelSync(const bool report) const
{
    if (!Pstream::parRun())
    {
        return false;
    }


    const polyBoundaryMesh& bm = *this;

    bool hasError = false;

    // Collect non-proc patches and check proc patches are last.
    wordList names(bm.size());
    wordList types(bm.size());

    label nonProci = 0;

    forAll(bm, patchi)
    {
        if (!isA<processorPolyPatch>(bm[patchi]))
        {
            if (nonProci != patchi)
            {
                // There is processor patch in between normal patches.
                hasError = true;

                if (debug || report)
                {
                    Pout<< " ***Problem with boundary patch " << patchi
                        << " named " << bm[patchi].name()
                        << " of type " <<  bm[patchi].type()
                        << ". The patch seems to be preceded by processor"
                        << " patches. This is can give problems."
                        << endl;
                }
            }
            else
            {
                names[nonProci] = bm[patchi].name();
                types[nonProci] = bm[patchi].type();
                nonProci++;
            }
        }
    }
    names.setSize(nonProci);
    types.setSize(nonProci);

    List<wordList> allNames(Pstream::nProcs());
    allNames[Pstream::myProcNo()] = names;
    Pstream::gatherList(allNames);
    Pstream::scatterList(allNames);

    List<wordList> allTypes(Pstream::nProcs());
    allTypes[Pstream::myProcNo()] = types;
    Pstream::gatherList(allTypes);
    Pstream::scatterList(allTypes);

    // Have every processor check but only master print error.

    for (label proci = 1; proci < allNames.size(); ++proci)
    {
        if
        (
            (allNames[proci] != allNames[0])
         || (allTypes[proci] != allTypes[0])
        )
        {
            hasError = true;

            if (debug || (report && Pstream::master()))
            {
                Info<< " ***Inconsistent patches across processors, "
                       "processor 0 has patch names:" << allNames[0]
                    << " patch types:" << allTypes[0]
                    << " processor " << proci << " has patch names:"
                    << allNames[proci]
                    << " patch types:" << allTypes[proci]
                    << endl;
            }
        }
    }

    return hasError;
}


bool Foam::polyBoundaryMesh::checkDefinition(const bool report) const
{
    const polyBoundaryMesh& bm = *this;

    bool hasError = false;

    HashSet<word> patchNames(2*size());

    forAll(bm, patchi)
    {
        if (!patchNames.insert(bm[patchi].name()) && !hasError)
        {
            hasError = true;

            Info<< " ****Duplicate boundary patch " << patchi
                << " named " << bm[patchi].name()
                << " of type " <<  bm[patchi].type()
                << "." << endl
                << "Suppressing future warnings." << endl;
        }
    }

    reduce(hasError, orOp<bool>());

    if (debug || report)
    {
        if (hasError)
        {
            Pout<< " ***Boundary definition is in error." << endl;
        }
        else
        {
            Info<< "    Boundary definition OK." << endl;
        }
    }

    return hasError;
}


void Foam::polyBoundaryMesh::movePoints(const pointField& p)
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchi)
        {
            operator[](patchi).initMovePoints(pBufs, p);
        }

        pBufs.finishedSends();

        forAll(*this, patchi)
        {
            operator[](patchi).movePoints(pBufs, p);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        forAll(patchSchedule, patchEvali)
        {
            const label patchi = patchSchedule[patchEvali].patch;

            if (patchSchedule[patchEvali].init)
            {
                operator[](patchi).initMovePoints(pBufs, p);
            }
            else
            {
                operator[](patchi).movePoints(pBufs, p);
            }
        }
    }
}


void Foam::polyBoundaryMesh::topoChange()
{
    nbrEdgesPtr_.clear();
    patchIndicesPtr_.clear();
    patchFaceIndicesPtr_.clear();
    groupPatchIndicesPtr_.clear();

    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchi)
        {
            operator[](patchi).initTopoChange(pBufs);
        }

        pBufs.finishedSends();

        forAll(*this, patchi)
        {
            operator[](patchi).topoChange(pBufs);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        forAll(patchSchedule, patchEvali)
        {
            const label patchi = patchSchedule[patchEvali].patch;

            if (patchSchedule[patchEvali].init)
            {
                operator[](patchi).initTopoChange(pBufs);
            }
            else
            {
                operator[](patchi).topoChange(pBufs);
            }
        }
    }
}


void Foam::polyBoundaryMesh::renamePatches
(
    const wordUList& newNames,
    const bool validBoundary
)
{
    polyPatchList& patches = *this;
    forAll(patches, patchi)
    {
        if (patches.set(patchi))
        {
            patches[patchi].rename(newNames);
        }
    }

    if (validBoundary)
    {
        topoChange();
    }
}


void Foam::polyBoundaryMesh::reorderPatches
(
    const labelUList& newToOld,
    const bool validBoundary
)
{
    // Change order of patches
    polyPatchList::shuffle(newToOld);

    // Adapt indices
    polyPatchList& patches = *this;
    forAll(patches, patchi)
    {
        if (patches.set(patchi))
        {
            patches[patchi].reorder(newToOld);
        }
    }

    if (validBoundary)
    {
        topoChange();
    }
}


bool Foam::polyBoundaryMesh::writeData(Ostream& os) const
{
    const polyPatchList& patches = *this;

    os  << patches.size() << nl << token::BEGIN_LIST << incrIndent << nl;

    forAll(patches, patchi)
    {
        os  << indent << patches[patchi].name() << nl
            << indent << token::BEGIN_BLOCK << nl
            << incrIndent << patches[patchi] << decrIndent
            << indent << token::END_BLOCK << endl;
    }

    os  << decrIndent << token::END_LIST;

    // Check state of IOstream
    os.check("polyBoundaryMesh::writeData(Ostream& os) const");

    return os.good();
}


bool Foam::polyBoundaryMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    return regIOobject::writeObject(fmt, ver, IOstream::UNCOMPRESSED, write);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

const Foam::polyPatch& Foam::polyBoundaryMesh::operator[]
(
    const word& patchName
) const
{
    const label patchi = findIndex(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch named " << patchName << " not found." << nl
            << "Available patch names: " << names() << endl
            << abort(FatalError);
    }

    return operator[](patchi);
}


Foam::polyPatch& Foam::polyBoundaryMesh::operator[]
(
    const word& patchName
)
{
    const label patchi = findIndex(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch named " << patchName << " not found." << nl
            << "Available patch names: " << names() << endl
            << abort(FatalError);
    }

    return operator[](patchi);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const polyBoundaryMesh& pbm)
{
    pbm.writeData(os);
    return os;
}


// ************************************************************************* //
