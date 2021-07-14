/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Cell layer addition/removal mesh modifier

\*---------------------------------------------------------------------------*/

#include "layerAdditionRemoval.H"
#include "polyTopoChanger.H"
#include "polyMesh.H"
#include "Time.H"
#include "primitiveMesh.H"
#include "polyTopoChange.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(layerAdditionRemoval, 0);
    addToRunTimeSelectionTable
    (
        polyMeshModifier,
        layerAdditionRemoval,
        dictionary
    );
}


const Foam::scalar Foam::layerAdditionRemoval::addDelta_ = 0.3;
const Foam::scalar Foam::layerAdditionRemoval::removeDelta_ = 0.1;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::layerAdditionRemoval::checkDefinition()
{
    if (!faceZoneID_.active())
    {
        FatalErrorInFunction
            << "Master face zone named " << faceZoneID_.name()
            << " cannot be found."
            << abort(FatalError);
    }

    if
    (
        minLayerThickness_ < vSmall
     || maxLayerThickness_ < minLayerThickness_
    )
    {
        FatalErrorInFunction
            << "Incorrect layer thickness definition."
            << abort(FatalError);
    }

    label nFaces = topoChanger().mesh().faceZones()[faceZoneID_.index()].size();

    reduce(nFaces, sumOp<label>());

    if (nFaces == 0)
    {
        FatalErrorInFunction
            << "Face extrusion zone contains no faces. "
            << "Please check your mesh definition."
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "Cell layer addition/removal object " << name() << " :" << nl
            << "    faceZoneID: " << faceZoneID_ << endl;
    }
}


Foam::scalar Foam::layerAdditionRemoval::readOldThickness
(
    const dictionary& dict
)
{
    return dict.lookupOrDefault("oldLayerThickness", -1.0);
}


void Foam::layerAdditionRemoval::clearAddressing() const
{
    if (pointsPairingPtr_)
    {
        if (debug)
        {
            Pout<< "layerAdditionRemoval::clearAddressing()" << nl
                << "    clearing pointsPairingPtr_" << endl;
        }

        deleteDemandDrivenData(pointsPairingPtr_);
    }

    if (facesPairingPtr_)
    {
        if (debug)
        {
            Pout<< "layerAdditionRemoval::clearAddressing()" << nl
                << "    clearing facesPairingPtr_" << endl;
        }

        deleteDemandDrivenData(facesPairingPtr_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::layerAdditionRemoval::layerAdditionRemoval
(
    const word& name,
    const label index,
    const polyTopoChanger& ptc,
    const word& zoneName,
    const scalar minThickness,
    const scalar maxThickness,
    const Switch thicknessFromVolume
)
:
    polyMeshModifier(name, index, ptc, true),
    faceZoneID_(zoneName, ptc.mesh().faceZones()),
    minLayerThickness_(minThickness),
    maxLayerThickness_(maxThickness),
    thicknessFromVolume_(thicknessFromVolume),
    oldLayerThickness_(-1.0),
    pointsPairingPtr_(nullptr),
    facesPairingPtr_(nullptr),
    triggerRemoval_(-1),
    triggerAddition_(-1)
{
    checkDefinition();
}


Foam::layerAdditionRemoval::layerAdditionRemoval
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyTopoChanger& ptc
)
:
    polyMeshModifier(name, index, ptc, Switch(dict.lookup("active"))),
    faceZoneID_(dict.lookup("faceZoneName"), ptc.mesh().faceZones()),
    minLayerThickness_(dict.lookup<scalar>("minLayerThickness")),
    maxLayerThickness_(dict.lookup<scalar>("maxLayerThickness")),
    thicknessFromVolume_
    (
        dict.lookupOrDefault<Switch>("thicknessFromVolume", true)
    ),
    oldLayerThickness_(readOldThickness(dict)),
    pointsPairingPtr_(nullptr),
    facesPairingPtr_(nullptr),
    triggerRemoval_(-1),
    triggerAddition_(-1)
{
    checkDefinition();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::layerAdditionRemoval::~layerAdditionRemoval()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::layerAdditionRemoval::changeTopology() const
{
    // Protect from multiple calculation in the same time-step
    if (triggerRemoval_ > -1 || triggerAddition_ > -1)
    {
        return true;
    }

    // Go through all the cells in the master layer and calculate
    // approximate layer thickness as the ratio of the cell volume and
    // face area in the face zone.
    // Layer addition:
    //     When the max thickness exceeds the threshold, trigger refinement.
    // Layer removal:
    //     When the min thickness falls below the threshold, trigger removal.

    const polyMesh& mesh = topoChanger().mesh();

    const faceZone& fz = mesh.faceZones()[faceZoneID_.index()];
    const labelList& mc = fz.masterCells();

    const scalarField& V = mesh.cellVolumes();
    const vectorField& S = mesh.faceAreas();

    if (min(V) < -vSmall)
    {
        FatalErrorInFunction
            << "negative cell volume. Error in mesh motion before "
            << "topological change.\n V: " << V
            << abort(FatalError);
    }

    scalar avgDelta = 0;
    scalar minDelta = great;
    scalar maxDelta = 0;
    label nDelta = 0;

    if (thicknessFromVolume_)
    {
        // Thickness calculated from cell volume/face area
        forAll(fz, facei)
        {
            scalar curDelta = V[mc[facei]]/mag(S[fz[facei]]);
            avgDelta += curDelta;
            minDelta = min(minDelta, curDelta);
            maxDelta = max(maxDelta, curDelta);
        }

        nDelta = fz.size();
    }
    else
    {
        // Thickness calculated from edges on layer
        const Map<label>& meshZonePointMap = fz().meshPointMap();

        // Edges with only one point on zone
        forAll(mc, facei)
        {
            const cell& cFaces = mesh.cells()[mc[facei]];
            const edgeList cellEdges(cFaces.edges(mesh.faces()));

            forAll(cellEdges, i)
            {
                const edge& e = cellEdges[i];

                if (meshZonePointMap.found(e[0]))
                {
                    if (!meshZonePointMap.found(e[1]))
                    {
                        scalar curDelta = e.mag(mesh.points());
                        avgDelta += curDelta;
                        nDelta++;
                        minDelta = min(minDelta, curDelta);
                        maxDelta = max(maxDelta, curDelta);
                    }
                }
                else
                {
                    if (meshZonePointMap.found(e[1]))
                    {
                        scalar curDelta = e.mag(mesh.points());
                        avgDelta += curDelta;
                        nDelta++;
                        minDelta = min(minDelta, curDelta);
                        maxDelta = max(maxDelta, curDelta);
                    }
                }
            }
        }
    }

    reduce(minDelta, minOp<scalar>());
    reduce(maxDelta, maxOp<scalar>());
    reduce(avgDelta, sumOp<scalar>());
    reduce(nDelta, sumOp<label>());

    avgDelta /= nDelta;

    if (debug)
    {
        Pout<< "bool layerAdditionRemoval::changeTopology() const "
            << " for object " << name() << " : " << nl
            << "Layer thickness: min: " << minDelta
            << " max: " << maxDelta << " avg: " << avgDelta
            << " old thickness: " << oldLayerThickness_ << nl
            << "Removal threshold: " << minLayerThickness_
            << " addition threshold: " << maxLayerThickness_ << endl;
    }

    bool topologicalChange = false;

    // If the thickness is decreasing and crosses the min thickness,
    // trigger removal
    if (oldLayerThickness_ < 0)
    {
        if (debug)
        {
            Pout<< "First step. No addition/removal" << endl;
        }

        // No topological changes allowed before first mesh motion
        oldLayerThickness_ = avgDelta;

        topologicalChange = false;
    }
    else if (avgDelta < oldLayerThickness_)
    {
        // Layers moving towards removal
        if (minDelta < minLayerThickness_)
        {
            // Check layer pairing
            if (setLayerPairing())
            {
                // A mesh layer detected.  Check that collapse is valid
                if (validCollapse())
                {
                    // At this point, info about moving the old mesh
                    // in a way to collapse the cells in the removed
                    // layer is available.  Not sure what to do with
                    // it.

                    if (debug)
                    {
                        Pout<< "bool layerAdditionRemoval::changeTopology() "
                            << " const for object " << name() << " : "
                            << "Triggering layer removal" << endl;
                    }

                    triggerRemoval_ = mesh.time().timeIndex();

                    // Old thickness looses meaning.
                    // Set it up to indicate layer removal
                    oldLayerThickness_ = great;

                    topologicalChange = true;
                }
                else
                {
                    // No removal, clear addressing
                    clearAddressing();
                }
            }
        }
        else
        {
            oldLayerThickness_ = avgDelta;
        }
    }
    else
    {
        // Layers moving towards addition
        if (maxDelta > maxLayerThickness_)
        {
            if (debug)
            {
                Pout<< "bool layerAdditionRemoval::changeTopology() const "
                    << " for object " << name() << " : "
                    << "Triggering layer addition" << endl;
            }

            triggerAddition_ = mesh.time().timeIndex();

            // Old thickness looses meaning.
            // Set it up to indicate layer removal
            oldLayerThickness_ = 0;

            topologicalChange = true;
        }
        else
        {
            oldLayerThickness_ = avgDelta;
        }
    }

    return topologicalChange;
}


void Foam::layerAdditionRemoval::setRefinement(polyTopoChange& ref) const
{
    // Insert the layer addition/removal instructions
    // into the topological change

    if (triggerRemoval_ == topoChanger().mesh().time().timeIndex())
    {
        removeCellLayer(ref);

        // Clear addressing.  This also resets the addition/removal data
        if (debug)
        {
            Pout<< "layerAdditionRemoval::setRefinement(polyTopoChange&) "
                << "for object " << name() << " : "
                << "Clearing addressing after layer removal" << endl;
        }

        triggerRemoval_ = -1;
        clearAddressing();
    }

    if (triggerAddition_ == topoChanger().mesh().time().timeIndex())
    {
        addCellLayer(ref);

        // Clear addressing.  This also resets the addition/removal data
        if (debug)
        {
            Pout<< "layerAdditionRemoval::setRefinement(polyTopoChange&) "
                << "for object " << name() << " : "
                << "Clearing addressing after layer addition" << endl;
        }

        triggerAddition_ = -1;
        clearAddressing();
    }
}


void Foam::layerAdditionRemoval::updateMesh(const mapPolyMesh&)
{
    if (debug)
    {
        Pout<< "layerAdditionRemoval::updateMesh(const mapPolyMesh&) "
            << "for object " << name() << " : "
            << "Clearing addressing on external request";

        if (pointsPairingPtr_ || facesPairingPtr_)
        {
            Pout<< "Pointers set." << endl;
        }
        else
        {
            Pout<< "Pointers not set." << endl;
        }
    }

    // Mesh has changed topologically.  Update local topological data
    faceZoneID_.update(topoChanger().mesh().faceZones());

    clearAddressing();
}


void Foam::layerAdditionRemoval::setMinLayerThickness(const scalar t) const
{
    if (t < vSmall || maxLayerThickness_ < t)
    {
        FatalErrorInFunction
            << "Incorrect layer thickness definition."
            << abort(FatalError);
    }

    minLayerThickness_ = t;
}


void Foam::layerAdditionRemoval::setMaxLayerThickness(const scalar t) const
{
    if (t < minLayerThickness_)
    {
        FatalErrorInFunction
            << "Incorrect layer thickness definition."
            << abort(FatalError);
    }

    maxLayerThickness_ = t;
}


void Foam::layerAdditionRemoval::write(Ostream& os) const
{
    os  << nl << type() << nl
        << name()<< nl
        << faceZoneID_ << nl
        << minLayerThickness_ << nl
        << oldLayerThickness_ << nl
        << maxLayerThickness_ << nl
        << thicknessFromVolume_ << endl;
}


void Foam::layerAdditionRemoval::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type()
        << token::END_STATEMENT << nl
        << "    faceZoneName " << faceZoneID_.name()
        << token::END_STATEMENT << nl
        << "    minLayerThickness " << minLayerThickness_
        << token::END_STATEMENT << nl
        << "    maxLayerThickness " << maxLayerThickness_
        << token::END_STATEMENT << nl
        << "    thicknessFromVolume " << thicknessFromVolume_
        << token::END_STATEMENT << nl
        << "    oldLayerThickness " << oldLayerThickness_
        << token::END_STATEMENT << nl
        << "    active " << active()
        << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// ************************************************************************* //
