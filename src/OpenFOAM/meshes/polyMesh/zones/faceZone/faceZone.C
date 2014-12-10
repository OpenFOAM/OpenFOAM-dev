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

#include "faceZone.H"
#include "addToRunTimeSelectionTable.H"
#include "faceZoneMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"
#include "mapPolyMesh.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceZone, 0);
    defineRunTimeSelectionTable(faceZone, dictionary);
    addToRunTimeSelectionTable(faceZone, faceZone, dictionary);
}

const char* const Foam::faceZone::labelsName = "faceLabels";

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceZone::calcFaceZonePatch() const
{
    if (debug)
    {
        Info<< "void faceZone::calcFaceZonePatch() const : "
            << "Calculating primitive patch"
            << endl;
    }

    if (patchPtr_)
    {
        FatalErrorIn
        (
            "void faceZone::calcFaceZonePatch() const"
        )   << "primitive face zone patch already calculated"
            << abort(FatalError);
    }

    patchPtr_ =
        new primitiveFacePatch
        (
            faceList(size()),
            zoneMesh().mesh().points()
        );

    primitiveFacePatch& patch = *patchPtr_;

    const faceList& f = zoneMesh().mesh().faces();

    const labelList& addr = *this;
    const boolList& flip = flipMap();

    forAll(addr, faceI)
    {
        if (flip[faceI])
        {
            patch[faceI] = f[addr[faceI]].reverseFace();
        }
        else
        {
            patch[faceI] = f[addr[faceI]];
        }
    }

    if (debug)
    {
        Info<< "void faceZone::calcFaceZonePatch() const : "
            << "Finished calculating primitive patch"
            << endl;
    }
}


void Foam::faceZone::calcCellLayers() const
{
    if (debug)
    {
        Info<< "void Foam::faceZone::calcCellLayers() const : "
            << "calculating master cells"
            << endl;
    }

    // It is an error to attempt to recalculate edgeCells
    // if the pointer is already set
    if (masterCellsPtr_ || slaveCellsPtr_)
    {
        FatalErrorIn("void faceZone::calcCellLayers() const")
            << "cell layers already calculated"
            << abort(FatalError);
    }
    else
    {
        // Go through all the faces in the master zone.  Choose the
        // master or slave cell based on the face flip

        const labelList& own = zoneMesh().mesh().faceOwner();
        const labelList& nei = zoneMesh().mesh().faceNeighbour();

        const labelList& mf = *this;

        const boolList& faceFlip = flipMap();

        masterCellsPtr_ = new labelList(mf.size());
        labelList& mc = *masterCellsPtr_;

        slaveCellsPtr_ = new labelList(mf.size());
        labelList& sc = *slaveCellsPtr_;

        forAll(mf, faceI)
        {
            label ownCellI = own[mf[faceI]];
            label neiCellI =
            (
                zoneMesh().mesh().isInternalFace(mf[faceI])
              ? nei[mf[faceI]]
              : -1
            );

            if (!faceFlip[faceI])
            {
                // Face is oriented correctly, no flip needed
                mc[faceI] = neiCellI;
                sc[faceI] = ownCellI;
            }
            else
            {
                mc[faceI] = ownCellI;
                sc[faceI] = neiCellI;
            }
        }
        //Info<< "masterCells: " << mc << endl;
        //Info<< "slaveCells: " << sc << endl;
    }
}


void Foam::faceZone::checkAddressing() const
{
    if (size() != flipMap_.size())
    {
        FatalErrorIn("void Foam::faceZone::checkAddressing() const")
            << "Different sizes of the addressing and flipMap arrays.  "
            << "Size of addressing: " << size()
            << " size of flip map: " << flipMap_.size()
            << abort(FatalError);
    }

    const labelList& mf = *this;

    // Note: nFaces, nCells might not be set yet on mesh so use owner size
    const label nFaces = zoneMesh().mesh().faceOwner().size();

    bool hasWarned = false;
    forAll(mf, i)
    {
        if (!hasWarned && (mf[i] < 0 || mf[i] >= nFaces))
        {
            WarningIn("void Foam::faceZone::checkAddressing() const")
                << "Illegal face index " << mf[i] << " outside range 0.."
                << nFaces-1 << endl;
            hasWarned = true;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::faceZone::faceZone
(
    const word& name,
    const labelUList& addr,
    const boolList& fm,
    const label index,
    const faceZoneMesh& zm
)
:
    zone(name, addr, index),
    flipMap_(fm),
    zoneMesh_(zm),
    patchPtr_(NULL),
    masterCellsPtr_(NULL),
    slaveCellsPtr_(NULL),
    mePtr_(NULL)
{
    checkAddressing();
}


Foam::faceZone::faceZone
(
    const word& name,
    const Xfer<labelList>& addr,
    const Xfer<boolList>& fm,
    const label index,
    const faceZoneMesh& zm
)
:
    zone(name, addr, index),
    flipMap_(fm),
    zoneMesh_(zm),
    patchPtr_(NULL),
    masterCellsPtr_(NULL),
    slaveCellsPtr_(NULL),
    mePtr_(NULL)
{
    checkAddressing();
}


Foam::faceZone::faceZone
(
    const word& name,
    const dictionary& dict,
    const label index,
    const faceZoneMesh& zm
)
:
    zone(name, dict, this->labelsName, index),
    flipMap_(dict.lookup("flipMap")),
    zoneMesh_(zm),
    patchPtr_(NULL),
    masterCellsPtr_(NULL),
    slaveCellsPtr_(NULL),
    mePtr_(NULL)
{
    checkAddressing();
}


Foam::faceZone::faceZone
(
    const faceZone& fz,
    const labelUList& addr,
    const boolList& fm,
    const label index,
    const faceZoneMesh& zm
)
:
    zone(fz, addr, index),
    flipMap_(fm),
    zoneMesh_(zm),
    patchPtr_(NULL),
    masterCellsPtr_(NULL),
    slaveCellsPtr_(NULL),
    mePtr_(NULL)
{
    checkAddressing();
}


Foam::faceZone::faceZone
(
    const faceZone& fz,
    const Xfer<labelList>& addr,
    const Xfer<boolList>& fm,
    const label index,
    const faceZoneMesh& zm
)
:
    zone(fz, addr, index),
    flipMap_(fm),
    zoneMesh_(zm),
    patchPtr_(NULL),
    masterCellsPtr_(NULL),
    slaveCellsPtr_(NULL),
    mePtr_(NULL)
{
    checkAddressing();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceZone::~faceZone()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::faceZoneMesh& Foam::faceZone::zoneMesh() const
{
    return zoneMesh_;
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
        //labelList faceCells(size());
        //
        //const labelList& own = zoneMesh().mesh().faceOwner();
        //
        //const labelList& faceLabels = *this;
        //
        //forAll(faceCells, faceI)
        //{
        //    faceCells[faceI] = own[faceLabels[faceI]];
        //}
        //
        //mePtr_ =
        //    new labelList
        //    (
        //        operator()().meshEdges
        //        (
        //            zoneMesh().mesh().edges(),
        //            zoneMesh().mesh().cellEdges(),
        //            faceCells
        //        )
        //    );

        mePtr_ =
            new labelList
            (
                operator()().meshEdges
                (
                    zoneMesh().mesh().edges(),
                    zoneMesh().mesh().pointEdges()
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


void Foam::faceZone::updateMesh(const mapPolyMesh& mpm)
{
    clearAddressing();

    labelList newAddressing(size());
    boolList newFlipMap(flipMap_.size());
    label nFaces = 0;

    const labelList& faceMap = mpm.reverseFaceMap();

    forAll(*this, i)
    {
        const label faceI = operator[](i);

        if (faceMap[faceI] >= 0)
        {
            newAddressing[nFaces] = faceMap[faceI];
            newFlipMap[nFaces] = flipMap_[i];       // Keep flip map.
            nFaces++;
        }
    }

    newAddressing.setSize(nFaces);
    newFlipMap.setSize(nFaces);

    transfer(newAddressing);
    flipMap_.transfer(newFlipMap);
}


bool Foam::faceZone::checkDefinition(const bool report) const
{
    return zone::checkDefinition(zoneMesh().mesh().faces().size(), report);
}


bool Foam::faceZone::checkParallelSync(const bool report) const
{
    const polyMesh& mesh = zoneMesh().mesh();
    const polyBoundaryMesh& bm = mesh.boundaryMesh();

    bool hasError = false;


    // Check that zone faces are synced
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        boolList neiZoneFace(mesh.nFaces()-mesh.nInternalFaces(), false);
        boolList neiZoneFlip(mesh.nFaces()-mesh.nInternalFaces(), false);
        forAll(*this, i)
        {
            const label faceI = operator[](i);

            if (!mesh.isInternalFace(faceI))
            {
                neiZoneFace[faceI-mesh.nInternalFaces()] = true;
                neiZoneFlip[faceI-mesh.nInternalFaces()] = flipMap()[i];
            }
        }
        boolList myZoneFace(neiZoneFace);
        syncTools::swapBoundaryFaceList(mesh, neiZoneFace);
        boolList myZoneFlip(neiZoneFlip);
        syncTools::swapBoundaryFaceList(mesh, neiZoneFlip);

        forAll(*this, i)
        {
            const label faceI = operator[](i);
            const label patchI = bm.whichPatch(faceI);

            if (patchI != -1 && bm[patchI].coupled())
            {
                const label bFaceI = faceI-mesh.nInternalFaces();

                // Check face in zone on both sides
                if (myZoneFace[bFaceI] != neiZoneFace[bFaceI])
                {
                    hasError = true;

                    if (report)
                    {
                        Pout<< " ***Problem with faceZone " << index()
                            << " named " << name()
                            << ". Face " << faceI
                            << " on coupled patch "
                            << bm[patchI].name()
                            << " is not consistent with its coupled neighbour."
                            << endl;
                    }
                    else
                    {
                        // w/o report - can stop checking now
                        break;
                    }
                }
                else if (myZoneFlip[bFaceI] == neiZoneFlip[bFaceI])
                {
                    // Flip state should be opposite.
                    hasError = true;

                    if (report)
                    {
                        Pout<< " ***Problem with faceZone " << index()
                            << " named " << name()
                            << ". Face " << faceI
                            << " on coupled patch "
                            << bm[patchI].name()
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
        patchPtr_->movePoints(p);
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

    writeEntry(this->labelsName, os);
    flipMap().writeEntry("flipMap", os);

    os  << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const faceZone& zn)
{
    zn.write(os);
    os.check("Ostream& operator<<(Ostream&, const faceZone&");
    return os;
}


// ************************************************************************* //
