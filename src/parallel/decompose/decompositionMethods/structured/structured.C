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

#include "structured.H"
#include "addToRunTimeSelectionTable.H"
#include "FaceCellWave.H"
#include "topoDistanceData.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionMethods
{
    defineTypeNameAndDebug(structured, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        structured,
        decomposer
    );

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        structured,
        distributor
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionMethods::structured::structured
(
    const dictionary& decompositionDict,
    const dictionary& methodDict
)
:
    decompositionMethod(decompositionDict)
{
    methodDict.lookup("patches") >> patches_;

    dictionary methodDict_(methodDict);
    methodDict_.set("numberOfSubdomains", nDomains());
    method_ = decompositionMethod::NewDecomposer(methodDict_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::decompositionMethods::structured::decompose
(
    const polyMesh& mesh,
    const pointField& cellCentres,
    const scalarField& cellWeights
)
{
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

    // Subset the layer of cells next to the patch
    fvMeshSubset subsetter(dynamic_cast<const fvMesh&>(mesh));
    subsetter.setLargeCellSubset(patchCells);
    const fvMesh& subMesh = subsetter.subMesh();
    pointField subCc(cellCentres, subsetter.cellMap());
    scalarField subWeights(cellWeights, subsetter.cellMap());

    // Decompose the layer of cells
    labelList subDecomp(method_().decompose(subMesh, subCc, subWeights));


    // Transfer to final decomposition
    labelList finalDecomp(cellCentres.size(), -1);
    forAll(subDecomp, i)
    {
        finalDecomp[subsetter.cellMap()[i]] = subDecomp[i];
    }

    // Field on cells and faces.
    List<topoDistanceData> cellData(mesh.nCells());
    List<topoDistanceData> faceData(mesh.nFaces());

    // Start of changes
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
            patchData[nFaces] = topoDistanceData(finalDecomp[fc[i]], 0);
            nFaces++;
        }
    }

    // Propagate information inwards
    FaceCellWave<topoDistanceData> deltaCalc
    (
        mesh,
        patchFaces,
        patchData,
        faceData,
        cellData,
        mesh.globalData().nTotalCells()+1
    );

    // And extract
    bool haveWarned = false;
    forAll(finalDecomp, celli)
    {
        if (!cellData[celli].valid(deltaCalc.data()))
        {
            if (!haveWarned)
            {
                WarningInFunction
                    << "Did not visit some cells, e.g. cell " << celli
                    << " at " << mesh.cellCentres()[celli] << endl
                    << "Assigning  these cells to domain 0." << endl;
                haveWarned = true;
            }
            finalDecomp[celli] = 0;
        }
        else
        {
            finalDecomp[celli] = cellData[celli].data();
        }
    }

    return finalDecomp;
}


Foam::labelList Foam::decompositionMethods::structured::decompose
(
    const labelListList& globalPointPoints,
    const pointField& points,
    const scalarField& pointWeights
)
{
    NotImplemented;

    return labelList::null();
}


// ************************************************************************* //
