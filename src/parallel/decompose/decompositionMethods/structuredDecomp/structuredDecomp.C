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

#include "structuredDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "FaceCellWave.H"
#include "topoDistanceData.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(structuredDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        structuredDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::structuredDecomp::structuredDecomp(const dictionary& decompositionDict)
:
    decompositionMethod(decompositionDict),
    methodDict_(decompositionDict_.subDict(typeName + "Coeffs")),
    patches_(methodDict_.lookup("patches"))
{
    methodDict_.set("numberOfSubdomains", nDomains());
    method_ = decompositionMethod::New(methodDict_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::structuredDecomp::parallelAware() const
{
    return method_().parallelAware();
}


Foam::labelList Foam::structuredDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& cc,
    const scalarField& cWeights
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
    pointField subCc(cc, subsetter.cellMap());
    scalarField subWeights(cWeights, subsetter.cellMap());

    // Decompose the layer of cells
    labelList subDecomp(method_().decompose(subMesh, subCc, subWeights));


    // Transfer to final decomposition
    labelList finalDecomp(cc.size(), -1);
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
    forAll(finalDecomp, cellI)
    {
        if (!cellData[cellI].valid(deltaCalc.data()))
        {
            if (!haveWarned)
            {
                WarningIn("structuredDecomp::decompose(..)")
                    << "Did not visit some cells, e.g. cell " << cellI
                    << " at " << mesh.cellCentres()[cellI] << endl
                    << "Assigning  these cells to domain 0." << endl;
                haveWarned = true;
            }
            finalDecomp[cellI] = 0;
        }
        else
        {
            finalDecomp[cellI] = cellData[cellI].data();
        }
    }

    return finalDecomp;
}


Foam::labelList Foam::structuredDecomp::decompose
(
    const labelListList& globalPointPoints,
    const pointField& points,
    const scalarField& pointWeights
)
{
    notImplemented
    (
        "structuredDecomp::decompose\n"
        "(\n"
        "    const labelListList&,\n"
        "    const pointField&,\n"
        "    const scalarField&\n"
        ")\n"
    );

    return labelList::null();
}


// ************************************************************************* //
