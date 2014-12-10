/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "structuredRenumber.H"
#include "addToRunTimeSelectionTable.H"
#include "topoDistanceData.H"
#include "fvMeshSubset.H"
#include "FaceCellWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(structuredRenumber, 0);

    addToRunTimeSelectionTable
    (
        renumberMethod,
        structuredRenumber,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::structuredRenumber::structuredRenumber
(
    const dictionary& renumberDict
)
:
    renumberMethod(renumberDict),
    methodDict_(renumberDict.subDict(typeName + "Coeffs")),
    patches_(methodDict_.lookup("patches")),
    //nLayers_(readLabel(methodDict_.lookup("nLayers"))),
    depthFirst_(methodDict_.lookup("depthFirst")),
    method_(renumberMethod::New(methodDict_)),
    reverse_(methodDict_.lookup("reverse"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::structuredRenumber::renumber
(
    const polyMesh& mesh,
    const pointField& points
) const
{
    if (points.size() != mesh.nCells())
    {
        FatalErrorIn
        (
            "structuredDecomp::renumber(const polyMesh&, const pointField&)"
        )   << "Number of points " << points.size()
            << " should equal the number of cells " << mesh.nCells()
            << exit(FatalError);
    }

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const labelHashSet patchIDs(pbm.patchSet(patches_));

    label nFaces = 0;
    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        nFaces += pbm[iter.key()].size();
    }


    // Extract a submesh.
    labelHashSet patchCells(2*nFaces);
    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        const labelUList& fc = pbm[iter.key()].faceCells();
        forAll(fc, i)
        {
            patchCells.insert(fc[i]);
        }
    }

    label nTotalSeeds = returnReduce(patchCells.size(), sumOp<label>());

    label nTotalCells = mesh.globalData().nTotalCells();
    const label nLayers = nTotalCells/nTotalSeeds;

    Info<< type() << " : seeding " << nTotalSeeds
        << " cells on " << nLayers << " layers" << nl
        << endl;


    // Avoid subsetMesh, FaceCellWave going through proc boundaries
    bool oldParRun = Pstream::parRun();
    Pstream::parRun() = false;


    // Work array. Used here to temporarily store the original-to-ordered
    // index. Later on used to store the ordered-to-original.
    labelList orderedToOld(points.size(), -1);

    // Subset the layer of cells next to the patch
    {
        fvMeshSubset subsetter(dynamic_cast<const fvMesh&>(mesh));
        subsetter.setLargeCellSubset(patchCells);
        const fvMesh& subMesh = subsetter.subMesh();

        pointField subPoints(points, subsetter.cellMap());

        // Decompose the layer of cells
        labelList subOrder(method_().renumber(subMesh, subPoints));

        labelList subOrigToOrdered(invert(subOrder.size(), subOrder));

        // Transfer to final decomposition
        forAll(subOrder, i)
        {
            orderedToOld[subsetter.cellMap()[i]] = subOrigToOrdered[i];
        }
    }


    // Walk out.
    labelList patchFaces(nFaces);
    List<topoDistanceData> patchData(nFaces);
    nFaces = 0;
    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        const polyPatch& pp = pbm[iter.key()];
        const labelUList& fc = pp.faceCells();
        forAll(fc, i)
        {
            patchFaces[nFaces] = pp.start()+i;
            patchData[nFaces] = topoDistanceData
            (
                orderedToOld[fc[i]],// passive data: order of originating face
                0                   // distance: layer
            );
            nFaces++;
        }
    }

    // Field on cells and faces.
    List<topoDistanceData> cellData(mesh.nCells());
    List<topoDistanceData> faceData(mesh.nFaces());

    // Propagate information inwards
    FaceCellWave<topoDistanceData> deltaCalc
    (
        mesh,
        patchFaces,
        patchData,
        faceData,
        cellData,
        nTotalCells+1
    );


    Pstream::parRun() = oldParRun;


    // And extract.
    // Note that distance is distance from face so starts at 1.
    bool haveWarned = false;
    forAll(orderedToOld, cellI)
    {
        if (!cellData[cellI].valid(deltaCalc.data()))
        {
            if (!haveWarned)
            {
                WarningIn("structuredDecomp::renumber(..)")
                    << "Did not visit some cells, e.g. cell " << cellI
                    << " at " << mesh.cellCentres()[cellI] << endl
                    << "Assigning these cells to domain 0." << endl;
                haveWarned = true;
            }
            orderedToOld[cellI] = 0;
        }
        else
        {
            label layerI = cellData[cellI].distance();
            if (depthFirst_)
            {
                orderedToOld[nLayers*cellData[cellI].data()+layerI] = cellI;
            }
            else
            {
                orderedToOld[cellData[cellI].data()+nLayers*layerI] = cellI;
            }
        }
    }

    // Return furthest away cell first
    if (reverse_)
    {
        reverse(orderedToOld);
    }

    return orderedToOld;
}


// ************************************************************************* //
