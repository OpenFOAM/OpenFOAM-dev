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

#include "faceZone.H"
#include "faceZoneList.H"
#include "polyMesh.H"
#include "polyTopoChangeMap.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    typedef Zone<faceZone, faceZoneList> faceZoneType;
    defineTemplateRunTimeSelectionTable(faceZoneType, dictionary);

    defineTypeNameAndDebug(faceZone, 0);
    addToRunTimeSelectionTable(faceZone, faceZone, dictionary);
}

const char* const Foam::faceZone::labelsName = "faceLabels";


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::faceZone::calcFaceZonePatch() const
{
    if (debug)
    {
        InfoInFunction << "Calculating primitive patch" << endl;
    }

    if (patchPtr_)
    {
        FatalErrorInFunction
            << "primitive face zone patch already calculated"
            << abort(FatalError);
    }

    patchPtr_ =
        new primitiveFacePatch
        (
            faceList(size()),
            zones().mesh().points()
        );

    primitiveFacePatch& patch = *patchPtr_;

    const faceList& f = zones().mesh().faces();

    const labelList& addr = *this;
    const boolList& flip = flipMap();

    forAll(addr, facei)
    {
        if (flip[facei])
        {
            patch[facei] = f[addr[facei]].reverseFace();
        }
        else
        {
            patch[facei] = f[addr[facei]];
        }
    }

    if (debug)
    {
        InfoInFunction << "Finished calculating primitive patch" << endl;
    }
}


void Foam::faceZone::calcCellLayers() const
{
    if (debug)
    {
        InfoInFunction << "Calculating master cells" << endl;
    }

    // It is an error to attempt to recalculate edgeCells
    // if the pointer is already set
    if (masterCellsPtr_ || slaveCellsPtr_)
    {
        FatalErrorInFunction
            << "cell layers already calculated"
            << abort(FatalError);
    }
    else
    {
        // Go through all the faces in the master zone.  Choose the
        // master or slave cell based on the face flip

        const labelList& own = zones().mesh().faceOwner();
        const labelList& nei = zones().mesh().faceNeighbour();

        const labelList& mf = *this;

        const boolList& faceFlip = flipMap();

        masterCellsPtr_ = new labelList(mf.size());
        labelList& mc = *masterCellsPtr_;

        slaveCellsPtr_ = new labelList(mf.size());
        labelList& sc = *slaveCellsPtr_;

        forAll(mf, facei)
        {
            label ownCelli = own[mf[facei]];
            label neiCelli =
            (
                zones().mesh().isInternalFace(mf[facei])
              ? nei[mf[facei]]
              : -1
            );

            if (!faceFlip[facei])
            {
                // Face is oriented correctly, no flip needed
                mc[facei] = neiCelli;
                sc[facei] = ownCelli;
            }
            else
            {
                mc[facei] = ownCelli;
                sc[facei] = neiCelli;
            }
        }
    }
}


void Foam::faceZone::checkAddressing() const
{
    if (size() != flipMap_.size())
    {
        FatalErrorInFunction
            << "Size of addressing: " << size()
            << " size of flip map: " << flipMap_.size()
            << abort(FatalError);
    }

    const labelList& mf = *this;

    // Note: nFaces, nCells might not be set yet on mesh so use owner size
    const label nFaces = zones().mesh().faceOwner().size();

    bool hasWarned = false;
    forAll(mf, i)
    {
        if (!hasWarned && (mf[i] < 0 || mf[i] >= nFaces))
        {
            WarningInFunction
                << "Illegal face index " << mf[i] << " outside range 0.."
                << nFaces-1 << endl;
            hasWarned = true;
        }
    }
}


void Foam::faceZone::reset(const Map<bool>& newIndices)
{
    clearAddressing();

    labelList& indices = *this;

    indices.setSize(newIndices.size());
    flipMap_.setSize(newIndices.size());

    const List<Map<bool>::const_iterator> sortedNewIndices(newIndices.sorted());

    forAll(sortedNewIndices, i)
    {
        indices[i] = sortedNewIndices[i].key();
        flipMap_[i] = sortedNewIndices[i]();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceZone::faceZone
(
    const word& name,
    const labelUList& addr,
    const boolList& fm,
    const faceZoneList& mz
)
:
    Zone<faceZone, faceZoneList>(name, addr, mz),
    flipMap_(fm),
    patchPtr_(nullptr),
    masterCellsPtr_(nullptr),
    slaveCellsPtr_(nullptr),
    mePtr_(nullptr)
{
    checkAddressing();
}


Foam::faceZone::faceZone
(
    const word& name,
    labelList&& addr,
    boolList&& fm,
    const faceZoneList& mz
)
:
    Zone<faceZone, faceZoneList>(name, move(addr), mz),
    flipMap_(move(fm)),
    patchPtr_(nullptr),
    masterCellsPtr_(nullptr),
    slaveCellsPtr_(nullptr),
    mePtr_(nullptr)
{
    checkAddressing();
}


Foam::faceZone::faceZone
(
    const word& name,
    const dictionary& dict,
    const faceZoneList& mz
)
:
    Zone<faceZone, faceZoneList>(name, dict, mz),
    flipMap_(dict.lookup("flipMap")),
    patchPtr_(nullptr),
    masterCellsPtr_(nullptr),
    slaveCellsPtr_(nullptr),
    mePtr_(nullptr)
{
    checkAddressing();
}


Foam::faceZone::faceZone
(
    const faceZone& fz,
    const labelUList& addr,
    const boolList& fm,
    const faceZoneList& mz
)
:
    Zone<faceZone, faceZoneList>(fz, addr, mz),
    flipMap_(fm),
    patchPtr_(nullptr),
    masterCellsPtr_(nullptr),
    slaveCellsPtr_(nullptr),
    mePtr_(nullptr)
{
    checkAddressing();
}


Foam::faceZone::faceZone
(
    const faceZone& fz,
    labelList&& addr,
    boolList&& fm,
    const faceZoneList& mz
)
:
    Zone<faceZone, faceZoneList>(fz, move(addr), mz),
    flipMap_(move(fm)),
    patchPtr_(nullptr),
    masterCellsPtr_(nullptr),
    slaveCellsPtr_(nullptr),
    mePtr_(nullptr)
{
    checkAddressing();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceZone::~faceZone()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::primitiveFacePatch& Foam::faceZone::patch() const
{
    if (!patchPtr_)
    {
        calcFaceZonePatch();
    }

    return *patchPtr_;
}


const Foam::labelList& Foam::faceZone::masterCells() const
{
    if (!masterCellsPtr_)
    {
        calcCellLayers();
    }

    return *masterCellsPtr_;
}


const Foam::labelList& Foam::faceZone::slaveCells() const
{
    if (!slaveCellsPtr_)
    {
        calcCellLayers();
    }

    return *slaveCellsPtr_;
}


const Foam::labelList& Foam::faceZone::meshEdges() const
{
    if (!mePtr_)
    {
        mePtr_ =
            new labelList
            (
                patch().meshEdges
                (
                    zones().mesh().edges(),
                    zones().mesh().pointEdges()
                )
            );
    }

    return *mePtr_;
}


void Foam::faceZone::clearAddressing()
{
    Zone::clearAddressing();

    deleteDemandDrivenData(patchPtr_);

    deleteDemandDrivenData(masterCellsPtr_);
    deleteDemandDrivenData(slaveCellsPtr_);

    deleteDemandDrivenData(mePtr_);
}


void Foam::faceZone::resetAddressing
(
    const labelUList& addr,
    const boolList& flipMap
)
{
    clearAddressing();
    labelList::operator=(addr);
    flipMap_ = flipMap;
}


bool Foam::faceZone::checkDefinition(const bool report) const
{
    return Zone::checkDefinition
    (
        zones().mesh().faces().size(),
        report
    );
}


bool Foam::faceZone::checkParallelSync(const bool report) const
{
    const polyMesh& mesh = zones().mesh();
    const polyBoundaryMesh& bm = mesh.boundaryMesh();

    bool hasError = false;


    // Check that zone faces are synced
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        boolList neiZoneFace(mesh.nFaces()-mesh.nInternalFaces(), false);
        boolList neiZoneFlip(mesh.nFaces()-mesh.nInternalFaces(), false);
        forAll(*this, i)
        {
            const label facei = operator[](i);

            if (!mesh.isInternalFace(facei))
            {
                neiZoneFace[facei-mesh.nInternalFaces()] = true;
                neiZoneFlip[facei-mesh.nInternalFaces()] = flipMap()[i];
            }
        }
        boolList myZoneFace(neiZoneFace);
        syncTools::swapBoundaryFaceList(mesh, neiZoneFace);
        boolList myZoneFlip(neiZoneFlip);
        syncTools::swapBoundaryFaceList(mesh, neiZoneFlip);

        forAll(*this, i)
        {
            const label facei = operator[](i);
            const label patchi = bm.whichPatch(facei);

            if (patchi != -1 && bm[patchi].coupled())
            {
                const label bFacei = facei-mesh.nInternalFaces();

                // Check face in zone on both sides
                if (myZoneFace[bFacei] != neiZoneFace[bFacei])
                {
                    hasError = true;

                    if (report)
                    {
                        Pout<< " ***Problem with faceZone " << name()
                            << ". Face " << facei
                            << " on coupled patch "
                            << bm[patchi].name()
                            << " is not consistent with its coupled neighbour."
                            << endl;
                    }
                    else
                    {
                        // w/o report - can stop checking now
                        break;
                    }
                }
                else if (myZoneFlip[bFacei] == neiZoneFlip[bFacei])
                {
                    // Flip state should be opposite.
                    hasError = true;

                    if (report)
                    {
                        Pout<< " ***Problem with faceZone " << name()
                            << ". Face " << facei
                            << " on coupled patch "
                            << bm[patchi].name()
                            << " does not have consistent flipMap"
                            << " across coupled faces."
                            << endl;
                    }
                    else
                    {
                        // w/o report - can stop checking now
                        break;
                    }
                }
            }
        }
    }

    return returnReduce(hasError, orOp<bool>());
}


void Foam::faceZone::insert(const Map<bool>& newIndices)
{
    Map<bool> indices(*this, flipMap_);
    indices.insert(newIndices);
    reset(indices);
}


void Foam::faceZone::swap(faceZone& fz)
{
    Zone::swap(fz);
    flipMap_.swap(fz.flipMap_);
}


void Foam::faceZone::topoChange(const polyTopoChangeMap& map)
{
    Map<bool> indices;
    const labelList& faceMap = map.faceMap();
    const labelList& reverseFaceMap = map.reverseFaceMap();
    const labelHashSet& flipFaceFlux = map.flipFaceFlux();

    forAll(faceMap, facei)
    {
        const label i = localIndex(faceMap[facei]);
        if (faceMap[facei] >= 0 && i != -1)
        {
            indices.insert
            (
                facei,
                flipFaceFlux.found(facei)
              ? !flipMap_[i]
              : flipMap_[i]
            );
        }
    }

    forAll(reverseFaceMap, facei)
    {
        const label i = localIndex(facei);
        if (reverseFaceMap[facei] >= 0 && i != -1)
        {
            indices.insert
            (
                reverseFaceMap[facei],
                flipFaceFlux.found(reverseFaceMap[facei])
              ? !flipMap_[i]
              : flipMap_[i]
            );
        }
    }

    reset(indices);
}


void Foam::faceZone::movePoints(const pointField& p)
{
    if (patchPtr_)
    {
        patchPtr_->clearGeom();
    }
}


void Foam::faceZone::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT << nl;

    writeEntry(os, this->labelsName, *this);
    writeEntry(os, "flipMap", flipMap());

    os  << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::faceZone::operator=(const faceZone& zn)
{
    Zone::operator=(zn);
    flipMap_ = zn.flipMap_;
}


void Foam::faceZone::operator=(faceZone&& zn)
{
    Zone::operator=(move(zn));
    flipMap_ = move(zn.flipMap_);
}


// ************************************************************************* //
