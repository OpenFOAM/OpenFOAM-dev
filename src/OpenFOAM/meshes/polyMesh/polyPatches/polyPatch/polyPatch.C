/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "polyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "SubField.H"
#include "entry.H"
#include "dictionary.H"
#include "pointPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyPatch, 0);

    int polyPatch::disallowGenericPolyPatch
    (
        debug::debugSwitch("disallowGenericPolyPatch", 0)
    );

    defineRunTimeSelectionTable(polyPatch, word);
    defineRunTimeSelectionTable(polyPatch, dictionary);

    addToRunTimeSelectionTable(polyPatch, polyPatch, word);
    addToRunTimeSelectionTable(polyPatch, polyPatch, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::polyPatch::movePoints(PstreamBuffers&, const pointField& p)
{
    primitivePatch::movePoints(p);
}

void Foam::polyPatch::updateMesh(PstreamBuffers&)
{
    primitivePatch::clearGeom();
    clearAddressing();
}


void Foam::polyPatch::clearGeom()
{
    primitivePatch::clearGeom();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyPatch::polyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    patchIdentifier(name, index),
    primitivePatch
    (
        faceSubList(bm.mesh().faces(), size, start),
        bm.mesh().points()
    ),
    start_(start),
    boundaryMesh_(bm),
    faceCellsPtr_(NULL),
    mePtr_(NULL)
{
    if
    (
        patchType != word::null
     && constraintType(patchType)
     && findIndex(inGroups(), patchType) == -1
    )
    {
        inGroups().append(patchType);
    }
}


Foam::polyPatch::polyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    patchIdentifier(name, dict, index),
    primitivePatch
    (
        faceSubList
        (
            bm.mesh().faces(),
            readLabel(dict.lookup("nFaces")),
            readLabel(dict.lookup("startFace"))
        ),
        bm.mesh().points()
    ),
    start_(readLabel(dict.lookup("startFace"))),
    boundaryMesh_(bm),
    faceCellsPtr_(NULL),
    mePtr_(NULL)
{
    if
    (
        patchType != word::null
     && constraintType(patchType)
     && findIndex(inGroups(), patchType) == -1
    )
    {
        inGroups().append(patchType);
    }
}


Foam::polyPatch::polyPatch
(
    const polyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    patchIdentifier(pp),
    primitivePatch
    (
        faceSubList
        (
            bm.mesh().faces(),
            pp.size(),
            pp.start()
        ),
        bm.mesh().points()
    ),
    start_(pp.start()),
    boundaryMesh_(bm),
    faceCellsPtr_(NULL),
    mePtr_(NULL)
{}


Foam::polyPatch::polyPatch
(
    const polyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    patchIdentifier(pp, index),
    primitivePatch
    (
        faceSubList
        (
            bm.mesh().faces(),
            newSize,
            newStart
        ),
        bm.mesh().points()
    ),
    start_(newStart),
    boundaryMesh_(bm),
    faceCellsPtr_(NULL),
    mePtr_(NULL)
{}


Foam::polyPatch::polyPatch
(
    const polyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    patchIdentifier(pp, index),
    primitivePatch
    (
        faceSubList
        (
            bm.mesh().faces(),
            mapAddressing.size(),
            newStart
        ),
        bm.mesh().points()
    ),
    start_(newStart),
    boundaryMesh_(bm),
    faceCellsPtr_(NULL),
    mePtr_(NULL)
{}


Foam::polyPatch::polyPatch(const polyPatch& p)
:
    patchIdentifier(p),
    primitivePatch(p),
    start_(p.start_),
    boundaryMesh_(p.boundaryMesh_),
    faceCellsPtr_(NULL),
    mePtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyPatch::~polyPatch()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::polyPatch::constraintType(const word& pt)
{
    return pointPatchField<scalar>::pointPatchConstructorTablePtr_->found(pt);
}


Foam::wordList Foam::polyPatch::constraintTypes()
{
    wordList cTypes(dictionaryConstructorTablePtr_->size());

    label i = 0;

    for
    (
        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->begin();
        cstrIter != dictionaryConstructorTablePtr_->end();
        ++cstrIter
    )
    {
        if (constraintType(cstrIter.key()))
        {
            cTypes[i++] = cstrIter.key();
        }
    }

    cTypes.setSize(i);

    return cTypes;
}


const Foam::polyBoundaryMesh& Foam::polyPatch::boundaryMesh() const
{
    return boundaryMesh_;
}


const Foam::vectorField::subField Foam::polyPatch::faceCentres() const
{
    return patchSlice(boundaryMesh().mesh().faceCentres());
}


const Foam::vectorField::subField Foam::polyPatch::faceAreas() const
{
    return patchSlice(boundaryMesh().mesh().faceAreas());
}


// Return the patch face neighbour cell centres
Foam::tmp<Foam::vectorField> Foam::polyPatch::faceCellCentres() const
{
    tmp<vectorField> tcc(new vectorField(size()));
    vectorField& cc = tcc();

    // get reference to global cell centres
    const vectorField& gcc = boundaryMesh_.mesh().cellCentres();

    const labelUList& faceCells = this->faceCells();

    forAll(faceCells, facei)
    {
        cc[facei] = gcc[faceCells[facei]];
    }

    return tcc;
}


const Foam::labelUList& Foam::polyPatch::faceCells() const
{
    if (!faceCellsPtr_)
    {
        faceCellsPtr_ = new labelList::subList
        (
            patchSlice(boundaryMesh().mesh().faceOwner())
        );
    }

    return *faceCellsPtr_;
}


const Foam::labelList& Foam::polyPatch::meshEdges() const
{
    if (!mePtr_)
    {
        mePtr_ =
            new labelList
            (
                primitivePatch::meshEdges
                (
                    boundaryMesh().mesh().edges(),
                    boundaryMesh().mesh().pointEdges()
                )
            );
    }

    return *mePtr_;
}


void Foam::polyPatch::clearAddressing()
{
    primitivePatch::clearTopology();
    primitivePatch::clearPatchMeshAddr();
    deleteDemandDrivenData(faceCellsPtr_);
    deleteDemandDrivenData(mePtr_);
}


void Foam::polyPatch::write(Ostream& os) const
{
    os.writeKeyword("type") << type() << token::END_STATEMENT << nl;
    patchIdentifier::write(os);
    os.writeKeyword("nFaces") << size() << token::END_STATEMENT << nl;
    os.writeKeyword("startFace") << start() << token::END_STATEMENT << nl;
}


void Foam::polyPatch::initOrder(PstreamBuffers&, const primitivePatch&) const
{}


bool Foam::polyPatch::order
(
    PstreamBuffers&,
    const primitivePatch&,
    labelList& faceMap,
    labelList& rotation
) const
{
    // Nothing changed.
    return false;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::polyPatch::operator=(const polyPatch& p)
{
    clearAddressing();

    patchIdentifier::operator=(p);
    primitivePatch::operator=(p);
    start_ = p.start_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const polyPatch& p)
{
    p.write(os);
    os.check("Ostream& operator<<(Ostream& os, const polyPatch& p");
    return os;
}


// ************************************************************************* //
