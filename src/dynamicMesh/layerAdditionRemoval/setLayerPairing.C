/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

Description
    Remove a layer of cells and prepare addressing data

\*---------------------------------------------------------------------------*/

#include "layerAdditionRemoval.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "polyTopoChange.H"
#include "oppositeFace.H"
#include "polyTopoChanger.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::layerAdditionRemoval::setLayerPairing() const
{
    // Note:
    // This is also the most complex part of the topological change.
    // Therefore it will be calculated here and stored as temporary
    // data until the actual topological change, after which it will
    // be cleared.

    // Algorithm for point collapse
    // 1)  Go through the master cell layer and for every face of
    //     the face zone find the opposite face in the master cell.
    //     Check the direction of the opposite face and adjust as
    //     necessary.  Check other faces to find an edge defining
    //     relative orientation of the two faces and adjust the face
    //     as necessary.  Once the face is adjusted, record the
    //     addressing between the master and slave vertex layer.

    const polyMesh& mesh = topoChanger().mesh();

    const labelList& mc =
        mesh.faceZones()[faceZoneID_.index()].masterCells();

    const labelList& mf = mesh.faceZones()[faceZoneID_.index()];

    const boolList& mfFlip =
        mesh.faceZones()[faceZoneID_.index()].flipMap();

    const faceList& faces = mesh.faces();
    const cellList& cells = mesh.cells();

    // Grab the local faces from the master zone
    const faceList& mlf =
        mesh.faceZones()[faceZoneID_.index()]().localFaces();

    const labelList& meshPoints =
        mesh.faceZones()[faceZoneID_.index()]().meshPoints();

    // Create a list of points to collapse for every point of
    // the master patch
    if (pointsPairingPtr_ || facesPairingPtr_)
    {
        FatalErrorInFunction
            << "Problem with layer pairing data"
            << abort(FatalError);
    }

    pointsPairingPtr_ = new labelList(meshPoints.size(), -1);
    labelList& ptc = *pointsPairingPtr_;

    facesPairingPtr_ = new labelList(mf.size(), -1);
    labelList& ftc = *facesPairingPtr_;

    if (debug > 1)
    {
        Pout<< "meshPoints: " << meshPoints << nl
            << "localPoints: "
            << mesh.faceZones()[faceZoneID_.index()]().localPoints()
            << endl;
    }

    // For all faces, create the mapping
    label nPointErrors = 0;
    label nFaceErrors = 0;

    forAll(mf, facei)
    {
        // Get the local master face
        face curLocalFace = mlf[facei];

        // Flip face based on flip index to recover original orientation
        if (mfFlip[facei])
        {
            curLocalFace.flip();
        }

        // Get the opposing face from the master cell
        oppositeFace lidFace =
            cells[mc[facei]].opposingFace(mf[facei], faces);

        if (!lidFace.found())
        {
            // This is not a valid layer; cannot continue
            nFaceErrors++;
            continue;
        }

        if (debug > 1)
        {
            Pout<< "curMasterFace: " << faces[mf[facei]] << nl
                << "cell shape: " << mesh.cellShapes()[mc[facei]] << nl
                << "curLocalFace: " << curLocalFace << nl
                << "lidFace: " << lidFace
                << " master index: " << lidFace.masterIndex()
                << " oppositeIndex: " << lidFace.oppositeIndex() << endl;
        }

        // Grab the opposite face for face collapse addressing
        ftc[facei] = lidFace.oppositeIndex();

        // Using the local face insert the points into the lid list
        forAll(curLocalFace, pointi)
        {
            const label clp = curLocalFace[pointi];

            if (ptc[clp] == -1)
            {
                // Point not mapped yet.  Insert the label
                ptc[clp] = lidFace[pointi];
            }
            else
            {
                // Point mapped from some other face.  Check the label
                if (ptc[clp] != lidFace[pointi])
                {
                    nPointErrors++;

                    if (debug > 1)
                    {
                        Pout<< "Topological error in cell layer pairing.  "
                            << "This mesh is either topologically incorrect "
                            << "or the master face layer is not defined "
                            << "consistently.  Please check the "
                            << "face zone flip map." << nl
                            << "First index: " << ptc[clp]
                            << " new index: " << lidFace[pointi] << endl;
                    }
                }
            }
        }
    }

    reduce(nPointErrors, sumOp<label>());
    reduce(nFaceErrors, sumOp<label>());

    if (nPointErrors > 0 || nFaceErrors > 0)
    {
        clearAddressing();

        return false;
    }
    else
    {
        // Valid layer
        return true;
    }
}


const Foam::labelList& Foam::layerAdditionRemoval::pointsPairing() const
{
    if (!pointsPairingPtr_)
    {
        FatalErrorInFunction
            << "Problem with layer pairing data for object " << name()
            << abort(FatalError);
    }

    return *pointsPairingPtr_;
}

const Foam::labelList& Foam::layerAdditionRemoval::facesPairing() const
{
    if (!facesPairingPtr_)
    {
        FatalErrorInFunction
            << "Problem with layer pairing data for object " << name()
            << abort(FatalError);
    }

    return *facesPairingPtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::layerAdditionRemoval::modifyMotionPoints
(
    pointField& motionPoints
) const
{
    if (debug)
    {
        Pout<< "void layerAdditionRemoval::modifyMotionPoints("
            << "pointField& motionPoints) const for object "
            << name() << " : ";
    }

    if (debug)
    {
        Pout<< "No motion point adjustment" << endl;
    }
}


// ************************************************************************* //
