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

\*---------------------------------------------------------------------------*/

#include "CFCFaceToCellStencil.H"
#include "syncTools.H"
#include "emptyPolyPatch.H"
#include "dummyTransform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::CFCFaceToCellStencil::calcFaceBoundaryData
(
    labelListList& neiGlobal
) const
{
    const polyBoundaryMesh& patches = mesh().boundaryMesh();
    const label nBnd = mesh().nFaces()-mesh().nInternalFaces();
    const labelList& own = mesh().faceOwner();

    neiGlobal.setSize(nBnd);

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        label facei = pp.start();

        if (pp.coupled())
        {
            // For coupled faces get the faces of the cell on the other side
            forAll(pp, i)
            {
                const labelList& cFaces = mesh().cells()[own[facei]];

                labelList& globFaces = neiGlobal[facei-mesh().nInternalFaces()];
                globFaces.setSize(cFaces.size()-1);
                label globI = 0;

                forAll(cFaces, j)
                {
                    if (cFaces[j] != facei)
                    {
                        globFaces[globI++] = globalNumbering().toGlobal
                        (
                            cFaces[j]
                        );
                    }
                }
                facei++;
            }
        }
        else if (isA<emptyPolyPatch>(pp))
        {}
        else
        {
            // Do nothing since face itself already in stencil
        }
    }

    syncTools::syncBoundaryFaceList
    (
        mesh(),
        neiGlobal,
        eqOp<labelList>(),
        dummyTransform()
    );
}


void Foam::CFCFaceToCellStencil::calcCellStencil
(
    labelListList& globalCellFaces
) const
{
    const label nBnd = mesh().nFaces()-mesh().nInternalFaces();
    const labelList& own = mesh().faceOwner();
    const labelList& nei = mesh().faceNeighbour();


    // Calculate faces of coupled neighbour (in global numbering)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelListList neiGlobal(nBnd);
    calcFaceBoundaryData(neiGlobal);



    // Non-empty boundary faces
    boolList validBFace(mesh().nFaces()-mesh().nInternalFaces(), true);

    const polyBoundaryMesh& patches = mesh().boundaryMesh();
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (isA<emptyPolyPatch>(pp))
        {
            label bFacei = pp.start()-mesh().nInternalFaces();
            forAll(pp, i)
            {
                validBFace[bFacei++] = false;
            }
        }
    }


    // Determine faces of cellCells in global numbering
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    DynamicList<label> allGlobalFaces(100);

    globalCellFaces.setSize(mesh().nCells());
    forAll(globalCellFaces, celli)
    {
        const cell& cFaces = mesh().cells()[celli];

        allGlobalFaces.clear();

        // My faces first
        forAll(cFaces, i)
        {
            label facei = cFaces[i];

            if
            (
                mesh().isInternalFace(facei)
             || validBFace[facei-mesh().nInternalFaces()]
            )
            {
                allGlobalFaces.append(globalNumbering().toGlobal(facei));
            }
        }

        // faces of neighbouring cells second
        forAll(cFaces, i)
        {
            label facei = cFaces[i];

            if (mesh().isInternalFace(facei))
            {
                label nbrCelli = own[facei];
                if (nbrCelli == celli)
                {
                    nbrCelli = nei[facei];
                }
                const cell& nbrFaces = mesh().cells()[nbrCelli];

                forAll(nbrFaces, j)
                {
                    label nbrFacei = nbrFaces[j];

                    if
                    (
                        mesh().isInternalFace(nbrFacei)
                     || validBFace[nbrFacei-mesh().nInternalFaces()]
                    )
                    {
                        label nbrGlobalI = globalNumbering().toGlobal(nbrFacei);

                        // Check if already there. Note:should use hashset?
                        if (findIndex(allGlobalFaces, nbrGlobalI) == -1)
                        {
                            allGlobalFaces.append(nbrGlobalI);
                        }
                    }
                }
            }
            else
            {
                const labelList& nbrGlobalFaces =
                    neiGlobal[facei-mesh().nInternalFaces()];

                forAll(nbrGlobalFaces, j)
                {
                    label nbrGlobalI = nbrGlobalFaces[j];

                    // Check if already there. Note:should use hashset?
                    if (findIndex(allGlobalFaces, nbrGlobalI) == -1)
                    {
                        allGlobalFaces.append(nbrGlobalI);
                    }
                }
            }
        }

        globalCellFaces[celli] = allGlobalFaces;
        // Pout<< "** cell:" << celli
        //    << " at:" << mesh().cellCentres()[celli]
        //    << endl;
        // const labelList& globalFaces = globalCellFaces[celli];
        // forAll(globalFaces, i)
        //{
        //    label facei = globalNumbering().toLocal(globalFaces[i]);
        //    Pout<< "    face:" << facei
        //        << " at:" << mesh().faceCentres()[facei]
        //        << endl;
        //}
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CFCFaceToCellStencil::CFCFaceToCellStencil(const polyMesh& mesh)
:
    faceToCellStencil(mesh)
{
    // Calculate per cell the (face) connected cells (in global numbering)
    calcCellStencil(*this);
}


// ************************************************************************* //
