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

#include "faceZone.H"
#include "addToRunTimeSelectionTable.H"
#include "meshFaceZones.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceZone, 0);
    defineRunTimeSelectionTable(faceZone, dictionary);
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
            meshZones().mesh().points()
        );

    primitiveFacePatch& patch = *patchPtr_;

    const faceList& f = meshZones().mesh().faces();

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

        const labelList& own = meshZones().mesh().faceOwner();
        const labelList& nei = meshZones().mesh().faceNeighbour();

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
                meshZones().mesh().isInternalFace(mf[facei])
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
    const label nFaces = meshZones().mesh().faceOwner().size();

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceZone::faceZone
(
    const word& name,
    const labelUList& addr,
    const boolList& fm,
    const label index,
    const meshFaceZones& mz
)
:
    zone(name, addr, index),
    flipMap_(fm),
    meshZones_(mz),
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
    const label index,
    const meshFaceZones& mz
)
:
    zone(name, move(addr), index),
    flipMap_(move(fm)),
    meshZones_(mz),
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
    const label index,
    const meshFaceZones& mz
)
:
    zone(name, dict, this->labelsName, index),
    flipMap_(dict.lookup("flipMap")),
    meshZones_(mz),
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
    const label index,
    const meshFaceZones& mz
)
:
    zone(fz, addr, index),
    flipMap_(fm),
    meshZones_(mz),
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
    const label index,
    const meshFaceZones& mz
)
:
    zone(fz, move(addr), index),
    flipMap_(move(fm)),
    meshZones_(mz),
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

const Foam::meshFaceZones& Foam::faceZone::meshZones() const
{
    return meshZones_;
}


Foam::label Foam::faceZone::whichFace(const label globalFaceID) const
{
    return zone::localID(globalFaceID);
}


const Foam::primitiveFacePatch& Foam::faceZone::operator()() const
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
                operator()().meshEdges
                (
                    meshZones().mesh().edges(),
                    meshZones().mesh().pointEdges()
                )
            );
    }

    return *mePtr_;
}


void Foam::faceZone::clearAddressing()
{
    zone::clearAddressing();

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


void Foam::faceZone::topoChange(const polyTopoChangeMap& map)
{
    clearAddressing();

    labelList newAddressing(size());
    boolList newFlipMap(flipMap_.size());
    label nFaces = 0;

    const labelList& faceMap = map.reverseFaceMap();

    forAll(*this, i)
    {
        const label facei = operator[](i);

        if (faceMap[facei] >= 0)
        {
            newAddressing[nFaces] = faceMap[facei];
            newFlipMap[nFaces] = flipMap_[i];       // Keep flip map.
            nFaces++;
        }
    }

    newAddressing.setSize(nFaces);
    newFlipMap.setSize(nFaces);

    transfer(newAddressing);
    flipMap_.transfer(newFlipMap);
}


void Foam::faceZone::mapMesh(const polyMeshMap&)
{
    NotImplemented;
}


bool Foam::faceZone::checkDefinition(const bool report) const
{
    return zone::checkDefinition(meshZones().mesh().faces().size(), report);
}


bool Foam::faceZone::checkParallelSync(const bool report) const
{
    const polyMesh& mesh = meshZones().mesh();
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
                        Pout<< " ***Problem with faceZone " << index()
                            << " named " << name()
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
                        Pout<< " ***Problem with faceZone " << index()
                            << " named " << name()
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


void Foam::faceZone::movePoints(const pointField& p)
{
    if (patchPtr_)
    {
        patchPtr_->clearGeom();
    }
}

void Foam::faceZone::write(Ostream& os) const
{
    os  << nl << name()
        << nl << static_cast<const labelList&>(*this)
        << nl << flipMap();
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
    clearAddressing();
    zone::operator=(zn);
    flipMap_ = zn.flipMap_;
}


void Foam::faceZone::operator=(faceZone&& zn)
{
    clearAddressing();
    zone::operator=(move(zn));
    flipMap_ = move(zn.flipMap_);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const faceZone& zn)
{
    zn.write(os);
    os.check("Ostream& operator<<(Ostream&, const faceZone&");
    return os;
}


// ************************************************************************* //
