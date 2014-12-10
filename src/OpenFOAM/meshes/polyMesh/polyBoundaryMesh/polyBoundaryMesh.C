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
            WarningIn
            (
                "polyBoundaryMesh::polyBoundaryMesh\n"
                "(\n"
                "    const IOobject&,\n"
                "    const polyMesh&\n"
                ")"
            )   << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
                << " does not support automatic rereading."
                << endl;
        }


        polyPatchList& patches = *this;

        // Read polyPatchList
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        patches.setSize(patchEntries.size());

        forAll(patches, patchI)
        {
            patches.set
            (
                patchI,
                polyPatch::New
                (
                    patchEntries[patchI].keyword(),
                    patchEntries[patchI].dict(),
                    patchI,
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
            WarningIn
            (
                "polyBoundaryMesh::polyBoundaryMesh\n"
                "(\n"
                "    const IOobject&,\n"
                "    const polyMesh&\n"
                "    const polyPatchList&\n"
                ")"
            )   << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
                << " does not support automatic rereading."
                << endl;
        }

        polyPatchList& patches = *this;

        // Read polyPatchList
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        patches.setSize(patchEntries.size());

        forAll(patches, patchI)
        {
            patches.set
            (
                patchI,
                polyPatch::New
                (
                    patchEntries[patchI].keyword(),
                    patchEntries[patchI].dict(),
                    patchI,
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
        forAll (patches, patchI)
        {
            patches.set(patchI, ppl[patchI].clone(*this).ptr());
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyBoundaryMesh::~polyBoundaryMesh()
{}


void Foam::polyBoundaryMesh::clearGeom()
{
    forAll(*this, patchI)
    {
        operator[](patchI).clearGeom();
    }
}


void Foam::polyBoundaryMesh::clearAddressing()
{
    neighbourEdgesPtr_.clear();
    patchIDPtr_.clear();
    groupPatchIDsPtr_.clear();

    forAll(*this, patchI)
    {
        operator[](patchI).clearAddressing();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyBoundaryMesh::calcGeometry()
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        Pstream::defaultCommsType == Pstream::blocking
     || Pstream::defaultCommsType == Pstream::nonBlocking
    )
    {
        forAll(*this, patchI)
        {
            operator[](patchI).initGeometry(pBufs);
        }

        pBufs.finishedSends();

        forAll(*this, patchI)
        {
            operator[](patchI).calcGeometry(pBufs);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        forAll(patchSchedule, patchEvali)
        {
            const label patchI = patchSchedule[patchEvali].patch;

            if (patchSchedule[patchEvali].init)
            {
                operator[](patchI).initGeometry(pBufs);
            }
            else
            {
                operator[](patchI).calcGeometry(pBufs);
            }
        }
    }
}


const Foam::List<Foam::labelPairList>&
Foam::polyBoundaryMesh::neighbourEdges() const
{
    if (Pstream::parRun())
    {
        WarningIn("polyBoundaryMesh::neighbourEdges() const")
            << "Neighbour edge addressing not correct across parallel"
            << " boundaries." << endl;
    }

    if (!neighbourEdgesPtr_.valid())
    {
        neighbourEdgesPtr_.reset(new List<labelPairList>(size()));
        List<labelPairList>& neighbourEdges = neighbourEdgesPtr_();

        // Initialize.
        label nEdgePairs = 0;
        forAll(*this, patchI)
        {
            const polyPatch& pp = operator[](patchI);

            neighbourEdges[patchI].setSize(pp.nEdges() - pp.nInternalEdges());

            forAll(neighbourEdges[patchI], i)
            {
                labelPair& edgeInfo = neighbourEdges[patchI][i];

                edgeInfo[0] = -1;
                edgeInfo[1] = -1;
            }

            nEdgePairs += pp.nEdges() - pp.nInternalEdges();
        }

        // From mesh edge (expressed as a point pair so as not to construct
        // point addressing) to patch + relative edge index.
        HashTable<labelPair, edge, Hash<edge> > pointsToEdge(nEdgePairs);

        forAll(*this, patchI)
        {
            const polyPatch& pp = operator[](patchI);

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

                HashTable<labelPair, edge, Hash<edge> >::iterator fnd =
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
                            patchI,
                            edgei - pp.nInternalEdges()
                        )
                    );
                }
                else
                {
                    // Second occurrence. Store.
                    const labelPair& edgeInfo = fnd();

                    neighbourEdges[patchI][edgei - pp.nInternalEdges()] =
                        edgeInfo;

                    neighbourEdges[edgeInfo[0]][edgeInfo[1]]
                         = labelPair(patchI, edgei - pp.nInternalEdges());

                    // Found all two occurrences of this edge so remove from
                    // hash to save space. Note that this will give lots of
                    // problems if the polyBoundaryMesh is multiply connected.
                    pointsToEdge.erase(meshEdge);
                }
            }
        }

        if (pointsToEdge.size())
        {
            FatalErrorIn("polyBoundaryMesh::neighbourEdges() const")
                << "Not all boundary edges of patches match up." << nl
                << "Is the outside of your mesh multiply connected?"
                << abort(FatalError);
        }

        forAll(*this, patchI)
        {
            const polyPatch& pp = operator[](patchI);

            const labelPairList& nbrEdges = neighbourEdges[patchI];

            forAll(nbrEdges, i)
            {
                const labelPair& edgeInfo = nbrEdges[i];

                if (edgeInfo[0] == -1 || edgeInfo[1] == -1)
                {
                    label edgeI = pp.nInternalEdges() + i;
                    const edge& e = pp.edges()[edgeI];

                    FatalErrorIn("polyBoundaryMesh::neighbourEdges() const")
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

    return neighbourEdgesPtr_();
}


const Foam::labelList& Foam::polyBoundaryMesh::patchID() const
{
    if (!patchIDPtr_.valid())
    {
        patchIDPtr_.reset
        (
            new labelList
            (
                mesh_.nFaces()
              - mesh_.nInternalFaces()
            )
        );
        labelList& patchID = patchIDPtr_();

        const polyBoundaryMesh& bm = *this;

        forAll(bm, patchI)
        {
            label bFaceI = bm[patchI].start() - mesh_.nInternalFaces();
            forAll(bm[patchI], i)
            {
                patchID[bFaceI++] = patchI;
            }
        }
    }
    return patchIDPtr_();
}


const Foam::HashTable<Foam::labelList, Foam::word>&
Foam::polyBoundaryMesh::groupPatchIDs() const
{
    if (!groupPatchIDsPtr_.valid())
    {
        groupPatchIDsPtr_.reset(new HashTable<labelList, word>(10));
        HashTable<labelList, word>& groupPatchIDs = groupPatchIDsPtr_();

        const polyBoundaryMesh& bm = *this;

        forAll(bm, patchI)
        {
            const wordList& groups = bm[patchI].inGroups();

            forAll(groups, i)
            {
                const word& name = groups[i];

                if (findPatchID(name) != -1)
                {
                    WarningIn("polyBoundaryMesh::groupPatchIDs() const")
                        << "Patch " << bm[patchI].name()
                        << " specifies a group " << name
                        << " which is also a patch name."
                        << " This might give problems later on." << endl;
                }


                HashTable<labelList, word>::iterator iter = groupPatchIDs.find
                (
                    name
                );

                if (iter != groupPatchIDs.end())
                {
                    iter().append(patchI);
                }
                else
                {
                    groupPatchIDs.insert(name, labelList(1, patchI));
                }
            }
        }
    }
    return groupPatchIDsPtr_();
}


void Foam::polyBoundaryMesh::setGroup
(
    const word& groupName,
    const labelList& patchIDs
)
{
    groupPatchIDsPtr_.clear();

    polyPatchList& patches = *this;

    boolList donePatch(patches.size(), false);

    // Add to specified patches
    forAll(patchIDs, i)
    {
        label patchI = patchIDs[i];
        polyPatch& pp = patches[patchI];

        if (!pp.inGroup(groupName))
        {
            pp.inGroups().append(groupName);
        }
        donePatch[patchI] = true;
    }

    // Remove from other patches
    forAll(patches, patchI)
    {
        if (!donePatch[patchI])
        {
            polyPatch& pp = patches[patchI];

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


Foam::wordList Foam::polyBoundaryMesh::names() const
{
    const polyPatchList& patches = *this;

    wordList t(patches.size());

    forAll(patches, patchI)
    {
        t[patchI] = patches[patchI].name();
    }

    return t;
}


Foam::wordList Foam::polyBoundaryMesh::types() const
{
    const polyPatchList& patches = *this;

    wordList t(patches.size());

    forAll(patches, patchI)
    {
        t[patchI] = patches[patchI].type();
    }

    return t;
}


Foam::wordList Foam::polyBoundaryMesh::physicalTypes() const
{
    const polyPatchList& patches = *this;

    wordList t(patches.size());

    forAll(patches, patchI)
    {
        t[patchI] = patches[patchI].physicalType();
    }

    return t;
}


Foam::labelList Foam::polyBoundaryMesh::findIndices
(
    const keyType& key,
    const bool usePatchGroups
) const
{
    DynamicList<label> indices;

    if (!key.empty())
    {
        if (key.isPattern())
        {
            indices = findStrings(key, this->names());

            if (usePatchGroups && groupPatchIDs().size())
            {
                labelHashSet indexSet(indices);

                const wordList allGroupNames = groupPatchIDs().toc();
                labelList groupIDs = findStrings(key, allGroupNames);
                forAll(groupIDs, i)
                {
                    const word& grpName = allGroupNames[groupIDs[i]];
                    const labelList& patchIDs = groupPatchIDs()[grpName];
                    forAll(patchIDs, j)
                    {
                        if (indexSet.insert(patchIDs[j]))
                        {
                            indices.append(patchIDs[j]);
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

            if (usePatchGroups && groupPatchIDs().size())
            {
                const HashTable<labelList, word>::const_iterator iter =
                    groupPatchIDs().find(key);

                if (iter != groupPatchIDs().end())
                {
                    labelHashSet indexSet(indices);

                    const labelList& patchIDs = iter();
                    forAll(patchIDs, j)
                    {
                        if (indexSet.insert(patchIDs[j]))
                        {
                            indices.append(patchIDs[j]);
                        }
                    }
                }
            }
        }
    }

    return indices;
}


Foam::label Foam::polyBoundaryMesh::findIndex(const keyType& key) const
{
    if (!key.empty())
    {
        if (key.isPattern())
        {
            labelList indices = this->findIndices(key);

            // return first element
            if (!indices.empty())
            {
                return indices[0];
            }
        }
        else
        {
            forAll(*this, i)
            {
                if (key == operator[](i).name())
                {
                    return i;
                }
            }
        }
    }

    // not found
    return -1;
}


Foam::label Foam::polyBoundaryMesh::findPatchID(const word& patchName) const
{
    const polyPatchList& patches = *this;

    forAll(patches, patchI)
    {
        if (patches[patchI].name() == patchName)
        {
            return patchI;
        }
    }

    // Patch not found
    if (debug)
    {
        Pout<< "label polyBoundaryMesh::findPatchID(const word&) const"
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
        FatalErrorIn
        (
            "polyBoundaryMesh::whichPatch(const label faceIndex) const"
        )   << "given label " << faceIndex
            << " greater than the number of geometric faces " << mesh().nFaces()
            << abort(FatalError);
    }


    forAll(*this, patchI)
    {
        const polyPatch& bp = operator[](patchI);

        if
        (
            faceIndex >= bp.start()
         && faceIndex < bp.start() + bp.size()
        )
        {
            return patchI;
        }
    }

    // If not in any of above, it is trouble!
    FatalErrorIn
    (
        "label polyBoundaryMesh::whichPatch(const label faceIndex) const"
    )   << "Cannot find face " << faceIndex << " in any of the patches "
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
        labelList patchIDs = findStrings(patchName, allPatchNames);

        forAll(patchIDs, j)
        {
            ids.insert(patchIDs[j]);
        }

        if (patchIDs.empty())
        {
            if (usePatchGroups)
            {
                const wordList allGroupNames = groupPatchIDs().toc();

                // Regard as group name
                labelList groupIDs = findStrings(patchName, allGroupNames);

                forAll(groupIDs, i)
                {
                    const word& name = allGroupNames[groupIDs[i]];
                    const labelList& extraPatchIDs = groupPatchIDs()[name];

                    forAll(extraPatchIDs, extraI)
                    {
                        ids.insert(extraPatchIDs[extraI]);
                    }
                }

                if (groupIDs.empty() && warnNotFound)
                {
                    WarningIn
                    (
                        "polyBoundaryMesh::patchSet"
                        "(const wordReList&, const bool, const bool) const"
                    )   << "Cannot find any patch or group names matching "
                        << patchName
                        << endl;
                }
            }
            else if (warnNotFound)
            {
                WarningIn
                (
                    "polyBoundaryMesh::patchSet"
                    "(const wordReList&, const bool, const bool) const"
                )   << "Cannot find any patch names matching " << patchName
                    << endl;
            }
        }
    }

    return ids;
}


void Foam::polyBoundaryMesh::matchGroups
(
    const labelUList& patchIDs,
    wordList& groups,
    labelHashSet& nonGroupPatches
) const
{
    // Current matched groups
    DynamicList<word> matchedGroups(1);

    // Current set of unmatched patches
    nonGroupPatches = labelHashSet(patchIDs);

    const HashTable<labelList, word>& groupPatchIDs = this->groupPatchIDs();
    for
    (
        HashTable<labelList,word>::const_iterator iter =
            groupPatchIDs.begin();
        iter != groupPatchIDs.end();
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

    label nonProcI = 0;

    forAll(bm, patchI)
    {
        if (!isA<processorPolyPatch>(bm[patchI]))
        {
            if (nonProcI != patchI)
            {
                // There is processor patch in between normal patches.
                hasError = true;

                if (debug || report)
                {
                    Pout<< " ***Problem with boundary patch " << patchI
                        << " named " << bm[patchI].name()
                        << " of type " <<  bm[patchI].type()
                        << ". The patch seems to be preceeded by processor"
                        << " patches. This is can give problems."
                        << endl;
                }
            }
            else
            {
                names[nonProcI] = bm[patchI].name();
                types[nonProcI] = bm[patchI].type();
                nonProcI++;
            }
        }
    }
    names.setSize(nonProcI);
    types.setSize(nonProcI);

    List<wordList> allNames(Pstream::nProcs());
    allNames[Pstream::myProcNo()] = names;
    Pstream::gatherList(allNames);
    Pstream::scatterList(allNames);

    List<wordList> allTypes(Pstream::nProcs());
    allTypes[Pstream::myProcNo()] = types;
    Pstream::gatherList(allTypes);
    Pstream::scatterList(allTypes);

    // Have every processor check but only master print error.

    for (label procI = 1; procI < allNames.size(); ++procI)
    {
        if
        (
            (allNames[procI] != allNames[0])
         || (allTypes[procI] != allTypes[0])
        )
        {
            hasError = true;

            if (debug || (report && Pstream::master()))
            {
                Info<< " ***Inconsistent patches across processors, "
                       "processor 0 has patch names:" << allNames[0]
                    << " patch types:" << allTypes[0]
                    << " processor " << procI << " has patch names:"
                    << allNames[procI]
                    << " patch types:" << allTypes[procI]
                    << endl;
            }
        }
    }

    return hasError;
}


bool Foam::polyBoundaryMesh::checkDefinition(const bool report) const
{
    label nextPatchStart = mesh().nInternalFaces();
    const polyBoundaryMesh& bm = *this;

    bool hasError = false;

    HashSet<word> patchNames(2*size());

    forAll(bm, patchI)
    {
        if (bm[patchI].start() != nextPatchStart && !hasError)
        {
            hasError = true;

            Info<< " ****Problem with boundary patch " << patchI
                << " named " << bm[patchI].name()
                << " of type " <<  bm[patchI].type()
                << ". The patch should start on face no " << nextPatchStart
                << " and the patch specifies " << bm[patchI].start()
                << "." << endl
                << "Possibly consecutive patches have this same problem."
                << " Suppressing future warnings." << endl;
        }

        if (!patchNames.insert(bm[patchI].name()) && !hasError)
        {
            hasError = true;

            Info<< " ****Duplicate boundary patch " << patchI
                << " named " << bm[patchI].name()
                << " of type " <<  bm[patchI].type()
                << "." << endl
                << "Suppressing future warnings." << endl;
        }

        nextPatchStart += bm[patchI].size();
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
        Pstream::defaultCommsType == Pstream::blocking
     || Pstream::defaultCommsType == Pstream::nonBlocking
    )
    {
        forAll(*this, patchI)
        {
            operator[](patchI).initMovePoints(pBufs, p);
        }

        pBufs.finishedSends();

        forAll(*this, patchI)
        {
            operator[](patchI).movePoints(pBufs, p);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        forAll(patchSchedule, patchEvali)
        {
            const label patchI = patchSchedule[patchEvali].patch;

            if (patchSchedule[patchEvali].init)
            {
                operator[](patchI).initMovePoints(pBufs, p);
            }
            else
            {
                operator[](patchI).movePoints(pBufs, p);
            }
        }
    }
}


void Foam::polyBoundaryMesh::updateMesh()
{
    neighbourEdgesPtr_.clear();
    patchIDPtr_.clear();
    groupPatchIDsPtr_.clear();

    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        Pstream::defaultCommsType == Pstream::blocking
     || Pstream::defaultCommsType == Pstream::nonBlocking
    )
    {
        forAll(*this, patchI)
        {
            operator[](patchI).initUpdateMesh(pBufs);
        }

        pBufs.finishedSends();

        forAll(*this, patchI)
        {
            operator[](patchI).updateMesh(pBufs);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        forAll(patchSchedule, patchEvali)
        {
            const label patchI = patchSchedule[patchEvali].patch;

            if (patchSchedule[patchEvali].init)
            {
                operator[](patchI).initUpdateMesh(pBufs);
            }
            else
            {
                operator[](patchI).updateMesh(pBufs);
            }
        }
    }
}


void Foam::polyBoundaryMesh::reorder
(
    const labelUList& oldToNew,
    const bool validBoundary
)
{
    // Change order of patches
    polyPatchList::reorder(oldToNew);

    // Adapt indices
    polyPatchList& patches = *this;

    forAll(patches, patchI)
    {
        patches[patchI].index() = patchI;
    }

    if (validBoundary)
    {
        updateMesh();
    }
}


bool Foam::polyBoundaryMesh::writeData(Ostream& os) const
{
    const polyPatchList& patches = *this;

    os  << patches.size() << nl << token::BEGIN_LIST << incrIndent << nl;

    forAll(patches, patchI)
    {
        os  << indent << patches[patchI].name() << nl
            << indent << token::BEGIN_BLOCK << nl
            << incrIndent << patches[patchI] << decrIndent
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
    IOstream::compressionType cmp
) const
{
    return regIOobject::writeObject(fmt, ver, IOstream::UNCOMPRESSED);
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

const Foam::polyPatch& Foam::polyBoundaryMesh::operator[]
(
    const word& patchName
) const
{
    const label patchI = findPatchID(patchName);

    if (patchI < 0)
    {
        FatalErrorIn
        (
            "polyBoundaryMesh::operator[](const word&) const"
        )   << "Patch named " << patchName << " not found." << nl
            << "Available patch names: " << names() << endl
            << abort(FatalError);
    }

    return operator[](patchI);
}


Foam::polyPatch& Foam::polyBoundaryMesh::operator[]
(
    const word& patchName
)
{
    const label patchI = findPatchID(patchName);

    if (patchI < 0)
    {
        FatalErrorIn
        (
            "polyBoundaryMesh::operator[](const word&)"
        )   << "Patch named " << patchName << " not found." << nl
            << "Available patch names: " << names() << endl
            << abort(FatalError);
    }

    return operator[](patchI);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const polyBoundaryMesh& pbm)
{
    pbm.writeData(os);
    return os;
}


// ************************************************************************* //
