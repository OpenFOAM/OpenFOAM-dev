/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
#include "OppositeFaceCellWave.H"

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
    methodDict_(renumberDict.optionalSubDict(typeName + "Coeffs")),
    patches_(methodDict_.lookup("patches")),
    nLayers_(methodDict_.lookupOrDefault<label>("nLayers", labelMax)),
    depthFirst_(methodDict_.lookup("depthFirst")),
    method_(renumberMethod::New(methodDict_)),
    reverse_(methodDict_.lookup("reverse"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::structuredRenumber::layerLess::operator()
(
    const label a,
    const label b
)
{
    const topoDistanceData& ta = distance_[a];
    const topoDistanceData& tb = distance_[b];

    int dummy;

    if (ta.valid(dummy))
    {
        if (tb.valid(dummy))
        {
            if (depthFirst_)
            {
                if (ta.data() < tb.data())
                {
                    // Sort column first
                    return true;
                }
                else if (ta.data() > tb.data())
                {
                    return false;
                }
                else
                {
                    // Same column. Sort according to layer
                    return ta.distance() < tb.distance();
                }
            }
            else
            {
                if (ta.distance() < tb.distance())
                {
                    return true;
                }
                else if (ta.distance() > tb.distance())
                {
                    return false;
                }
                else
                {
                    // Same layer; sort according to current values
                    return ta.data() < tb.data();
                }
            }
        }
        else
        {
            return true;
        }
    }
    else
    {
        if (tb.valid(dummy))
        {
            return false;
        }
        else
        {
            // Both not valid; fall back to cell indices for sorting
            return order_[a] < order_[b];
        }
    }
}


Foam::labelList Foam::structuredRenumber::renumber
(
    const polyMesh& mesh,
    const pointField& points
) const
{
    if (points.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Number of points " << points.size()
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
        << " cells on (estimated) " << nLayers << " layers" << nl
        << endl;


    // Work array. Used here to temporarily store the original-to-ordered
    // index. Later on used to store the ordered-to-original.
    labelList orderedToOld(mesh.nCells(), -1);

    // Subset the layer of cells next to the patch
    {
        fvMeshSubset subsetter(dynamic_cast<const fvMesh&>(mesh));
        subsetter.setLargeCellSubset(patchCells);
        const fvMesh& subMesh = subsetter.subMesh();

        pointField subPoints(points, subsetter.cellMap());

        // Locally renumber the layer of cells
        labelList subOrder(method_().renumber(subMesh, subPoints));

        labelList subOrigToOrdered(invert(subOrder.size(), subOrder));

        globalIndex globalSubCells(subOrder.size());

        // Transfer to final decomposition and convert into global numbering
        forAll(subOrder, i)
        {
            orderedToOld[subsetter.cellMap()[i]] =
                globalSubCells.toGlobal(subOrigToOrdered[i]);
        }
    }


    // Walk sub-ordering (=column index) out.
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
                orderedToOld[fc[i]],// passive data: global column
                0                   // distance: layer
            );
            nFaces++;
        }
    }

    // Field on cells and faces.
    List<topoDistanceData> cellData(mesh.nCells());
    List<topoDistanceData> faceData(mesh.nFaces());

    // Propagate information inwards
    OppositeFaceCellWave<topoDistanceData> deltaCalc
    (
        mesh,
        patchFaces,
        patchData,
        faceData,
        cellData,
        0
    );

    deltaCalc.iterate(nLayers_);

    Info<< type() << " : did not visit "
        << deltaCalc.getUnsetCells()
        << " cells out of " << nTotalCells
        << "; using " << method_().type() << " renumbering for these" << endl;

    // Get cell order using the method(). These values will get overwritten
    // by any visited cell so are used only if the number of nLayers is limited.
    labelList oldToOrdered
    (
        invert
        (
            mesh.nCells(),
            method_().renumber(mesh, points)
        )
    );

    // Use specialised sorting to sorted either layers or columns first
    // Done so that at no point we need to combine both into a single
    // index and we might run out of label size.
    sortedOrder
    (
        cellData,
        orderedToOld,
        layerLess(depthFirst_, oldToOrdered, cellData)
    );

    // Return furthest away cell first
    if (reverse_)
    {
        reverse(orderedToOld);
    }

    return orderedToOld;
}


// ************************************************************************* //
