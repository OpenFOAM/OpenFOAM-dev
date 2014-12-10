/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "refinementIterator.H"
#include "polyMesh.H"
#include "Time.H"
#include "refineCell.H"
#include "undoableMeshCutter.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "cellCuts.H"
#include "OFstream.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(refinementIterator, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::refinementIterator::refinementIterator
(
    polyMesh& mesh,
    undoableMeshCutter& meshRefiner,
    const cellLooper& cellWalker,
    const bool writeMesh
)
:
    edgeVertex(mesh),
    mesh_(mesh),
    meshRefiner_(meshRefiner),
    cellWalker_(cellWalker),
    writeMesh_(writeMesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::refinementIterator::~refinementIterator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Map<Foam::label> Foam::refinementIterator::setRefinement
(
    const List<refineCell>& refCells
)
{
    Map<label> addedCells(2*refCells.size());

    Time& runTime = const_cast<Time&>(mesh_.time());

    label nRefCells = refCells.size();

    label oldRefCells = -1;

    // Operate on copy.
    List<refineCell> currentRefCells(refCells);

    bool stop = false;

    do
    {
        if (writeMesh_)
        {
            // Need different times to write meshes.
            runTime++;
        }

        polyTopoChange meshMod(mesh_);

        if (debug)
        {
            Pout<< "refinementIterator : refining "
                << currentRefCells.size() << " cells." << endl;
        }

        // Determine cut pattern.
        cellCuts cuts(mesh_, cellWalker_, currentRefCells);

        label nCuts = cuts.nLoops();
        reduce(nCuts, sumOp<label>());

        if (nCuts == 0)
        {
            if (debug)
            {
                Pout<< "refinementIterator : exiting iteration since no valid"
                    << " loops found for " << currentRefCells.size()
                    << " cells" << endl;


                fileName cutsFile("failedCuts_" + runTime.timeName() + ".obj");

                Pout<< "Writing cuts for time " <<  runTime.timeName()
                    << " to " << cutsFile << endl;

                OFstream cutsStream(cutsFile);


                labelList refCells(currentRefCells.size());
                forAll(currentRefCells, i)
                {
                    refCells[i] = currentRefCells[i].cellNo();
                }
                meshTools::writeOBJ
                (
                    cutsStream,
                    mesh().cells(),
                    mesh().faces(),
                    mesh().points(),
                    refCells
                );
            }

            break;
        }

        if (debug)
        {
            fileName cutsFile("cuts_" + runTime.timeName() + ".obj");

            Pout<< "Writing cuts for time " <<  runTime.timeName()
                << " to " << cutsFile << endl;

            OFstream cutsStream(cutsFile);
            cuts.writeOBJ(cutsStream);
        }


        // Insert mesh refinement into polyTopoChange.
        meshRefiner_.setRefinement(cuts, meshMod);


        //
        // Do all changes
        //

        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh
        (
            mesh_,
            false
        );

        // Move mesh (since morphing does not do this)
        if (morphMap().hasMotionPoints())
        {
            mesh_.movePoints(morphMap().preMotionPoints());
        }

        // Update stored refinement pattern
        meshRefiner_.updateMesh(morphMap());

        // Write resulting mesh
        if (writeMesh_)
        {
            if (debug)
            {
                Pout<< "Writing refined polyMesh to time "
                    << runTime.timeName() << endl;
            }

            mesh_.write();
        }

        // Update currentRefCells for new cell numbers. Use helper function
        // in meshCutter class.
        updateLabels
        (
            morphMap->reverseCellMap(),
            currentRefCells
        );

        // Update addedCells for new cell numbers
        updateLabels
        (
            morphMap->reverseCellMap(),
            addedCells
        );

        // Get all added cells from cellCutter (already in new numbering
        // from meshRefiner.updateMesh call) and add to global list of added
        const Map<label>& addedNow = meshRefiner_.addedCells();

        forAllConstIter(Map<label>, addedNow, iter)
        {
            if (!addedCells.insert(iter.key(), iter()))
            {
                FatalErrorIn("refinementIterator")
                    << "Master cell " << iter.key()
                    << " already has been refined" << endl
                    << "Added cell:" << iter() << abort(FatalError);
            }
        }


        // Get failed refinement in new cell numbering and reconstruct input
        // to the meshRefiner. Is done by removing all refined cells from
        // current list of cells to refine.

        // Update refCells for new cell numbers.
        updateLabels
        (
            morphMap->reverseCellMap(),
            currentRefCells
        );

        // Pack refCells acc. to refined status
        nRefCells = 0;

        forAll(currentRefCells, refI)
        {
            const refineCell& refCell = currentRefCells[refI];

            if (!addedNow.found(refCell.cellNo()))
            {
                if (nRefCells != refI)
                {
                    currentRefCells[nRefCells++] =
                        refineCell
                        (
                            refCell.cellNo(),
                            refCell.direction()
                        );
                }
            }
        }

        oldRefCells = currentRefCells.size();

        currentRefCells.setSize(nRefCells);

        if (debug)
        {
            Pout<< endl;
        }

        // Stop only if all finished or all can't refine any further.
        stop = (nRefCells == 0) || (nRefCells == oldRefCells);
        reduce(stop, andOp<bool>());
    }
    while (!stop);


    if (nRefCells == oldRefCells)
    {
        WarningIn("refinementIterator")
            << "stopped refining."
            << "Did not manage to refine a single cell" << endl
            << "Wanted :" << oldRefCells << endl;
    }

    return addedCells;
}



// ************************************************************************* //
