/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "MGridGenGAMGAgglomeration.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "processorLduInterface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MGridGenGAMGAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        MGridGenGAMGAgglomeration,
        lduMesh
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::MGridGenGAMGAgglomeration::swap
(
    const lduInterfacePtrsList& interfaces,
    const labelUList& cellValues,
    PtrList<labelList>& nbrValues
) const
{
    // Initialise transfer of restrict addressing on the interface
    forAll(interfaces, inti)
    {
        if (interfaces.set(inti))
        {
            interfaces[inti].initInternalFieldTransfer
            (
                Pstream::commsTypes::nonBlocking,
                cellValues
            );
        }
    }

    if (Pstream::parRun())
    {
        Pstream::waitRequests();
    }

    // Get the interface agglomeration
    nbrValues.setSize(interfaces.size());
    forAll(interfaces, inti)
    {
        if (interfaces.set(inti))
        {
            nbrValues.set
            (
                inti,
                new labelList
                (
                    interfaces[inti].internalFieldTransfer
                    (
                        Pstream::commsTypes::nonBlocking,
                        cellValues
                    )
                )
            );
        }
    }
}


void Foam::MGridGenGAMGAgglomeration::getNbrAgglom
(
    const lduAddressing& addr,
    const lduInterfacePtrsList& interfaces,
    const PtrList<labelList>& nbrGlobalAgglom,
    labelList& cellToNbrAgglom
) const
{
    cellToNbrAgglom.setSize(addr.size());
    cellToNbrAgglom = -1;

    forAll(interfaces, inti)
    {
        if (interfaces.set(inti))
        {
            if (isA<processorLduInterface>(interfaces[inti]))
            {
                const processorLduInterface& pldui =
                    refCast<const processorLduInterface>(interfaces[inti]);

                if (pldui.myProcNo() > pldui.neighbProcNo())
                {
                    const labelUList& faceCells =
                        interfaces[inti].faceCells();
                    const labelList& nbrData = nbrGlobalAgglom[inti];

                    forAll(faceCells, i)
                    {
                        cellToNbrAgglom[faceCells[i]] = nbrData[i];
                    }
                }
            }
        }
    }
}


void Foam::MGridGenGAMGAgglomeration::detectSharedFaces
(
    const lduMesh& mesh,
    const labelList& value,
    labelHashSet& sharedFaces
) const
{
    const lduAddressing& addr = mesh.lduAddr();
    const labelUList& lower = addr.lowerAddr();
    const labelUList& upper = addr.upperAddr();

    sharedFaces.clear();
    sharedFaces.resize(addr.lowerAddr().size()/100);

    // Detect any faces in between same value
    forAll(lower, facei)
    {
        label lowerData = value[lower[facei]];
        label upperData = value[upper[facei]];

        if (lowerData != -1 && lowerData == upperData)
        {
            sharedFaces.insert(facei);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MGridGenGAMGAgglomeration::MGridGenGAMGAgglomeration
(
    const lduMesh& mesh,
    const dictionary& controlDict
)
:
    GAMGAgglomeration(mesh, controlDict),
    fvMesh_(refCast<const fvMesh>(mesh))
{
    // Min, max size of agglomerated cells
    label minSize(readLabel(controlDict.lookup("minSize")));
    label maxSize(readLabel(controlDict.lookup("maxSize")));

    // Number of iterations applied to improve agglomeration consistency across
    // processor boundaries
    label nProcConsistencyIter
    (
        readLabel(controlDict.lookup("nProcConsistencyIter"))
    );

    // Start geometric agglomeration from the cell volumes and areas of the mesh
    scalarField* VPtr = const_cast<scalarField*>(&fvMesh_.cellVolumes());

    scalarField magFaceAreas(sqrt(3.0)*mag(fvMesh_.faceAreas()));
    SubField<scalar> magSf(magFaceAreas, fvMesh_.nInternalFaces());

    scalarField* magSfPtr = const_cast<scalarField*>
    (
        &magSf.operator const scalarField&()
    );

    // Create the boundary area cell field
    scalarField* magSbPtr(new scalarField(fvMesh_.nCells(), 0));

    {
        scalarField& magSb = *magSbPtr;

        const labelList& own = fvMesh_.faceOwner();
        const vectorField& Sf = fvMesh_.faceAreas();

        forAll(Sf, facei)
        {
            if (!fvMesh_.isInternalFace(facei))
            {
                magSb[own[facei]] += mag(Sf[facei]);
            }
        }
    }

    // Agglomerate until the required number of cells in the coarsest level
    // is reached

    label nCreatedLevels = 0;

    while (nCreatedLevels < maxLevels_ - 1)
    {
        label nCoarseCells = -1;

        tmp<labelField> finalAgglomPtr = agglomerate
        (
            nCoarseCells,
            minSize,
            maxSize,
            meshLevel(nCreatedLevels).lduAddr(),
            *VPtr,
            *magSfPtr,
            *magSbPtr
        );

        // Adjust weights only
        for (int i=0; i<nProcConsistencyIter; i++)
        {
            const lduMesh& mesh = meshLevel(nCreatedLevels);
            const lduAddressing& addr = mesh.lduAddr();
            const lduInterfacePtrsList interfaces = mesh.interfaces();

            const labelField& agglom = finalAgglomPtr();

            // Global nubmering
            const globalIndex globalNumbering(nCoarseCells);

            labelField globalAgglom(addr.size());
            forAll(agglom, celli)
            {
                globalAgglom[celli] = globalNumbering.toGlobal(agglom[celli]);
            }

            // Get the interface agglomeration
            PtrList<labelList> nbrGlobalAgglom;
            swap(interfaces, globalAgglom, nbrGlobalAgglom);


            // Get the interface agglomeration on a cell basis (-1 for all
            // other cells)
            labelList cellToNbrAgglom;
            getNbrAgglom(addr, interfaces, nbrGlobalAgglom, cellToNbrAgglom);


            // Mark all faces in between cells with same nbragglomeration
            labelHashSet sharedFaces(addr.size()/100);
            detectSharedFaces(mesh, cellToNbrAgglom, sharedFaces);


            //- Note: in-place update of weights is more effective it seems?
            //        Should not be. fluke?
            //scalarField weights(*faceWeightsPtr);
            scalarField weights = *magSfPtr;
            forAllConstIter(labelHashSet, sharedFaces, iter)
            {
                label facei= iter.key();
                weights[facei] *= 2.0;
            }

            // Redo the agglomeration using the new weights
            finalAgglomPtr = agglomerate
            (
                nCoarseCells,
                minSize,
                maxSize,
                meshLevel(nCreatedLevels).lduAddr(),
                *VPtr,
                weights,
                *magSbPtr
            );
        }

        if (continueAgglomerating(finalAgglomPtr().size(), nCoarseCells))
        {
            nCells_[nCreatedLevels] = nCoarseCells;
            restrictAddressing_.set(nCreatedLevels, finalAgglomPtr);
        }
        else
        {
            break;
        }

        agglomerateLduAddressing(nCreatedLevels);

        // Agglomerate the cell volumes field for the next level
        {
            scalarField* aggVPtr
            (
                new scalarField(meshLevels_[nCreatedLevels].size())
            );

            // Restrict but no parallel agglomeration (not supported)
            restrictField(*aggVPtr, *VPtr, nCreatedLevels, false);

            if (nCreatedLevels)
            {
                delete VPtr;
            }

            VPtr = aggVPtr;
        }

        // Agglomerate the face areas field for the next level
        {
            scalarField* aggMagSfPtr
            (
                new scalarField
                (
                    meshLevels_[nCreatedLevels].upperAddr().size(),
                    0
                )
            );

            restrictFaceField(*aggMagSfPtr, *magSfPtr, nCreatedLevels);

            if (nCreatedLevels)
            {
                delete magSfPtr;
            }

            magSfPtr = aggMagSfPtr;
        }

        // Agglomerate the cell boundary areas field for the next level
        {
            scalarField* aggMagSbPtr
            (
                new scalarField(meshLevels_[nCreatedLevels].size())
            );

            // Restrict but no parallel agglomeration (not supported)
            restrictField(*aggMagSbPtr, *magSbPtr, nCreatedLevels, false);

            delete magSbPtr;
            magSbPtr = aggMagSbPtr;
        }

        nCreatedLevels++;
    }

    // Shrink the storage of the levels to those created
    compactLevels(nCreatedLevels);

    // Delete temporary geometry storage
    if (nCreatedLevels)
    {
        delete VPtr;
        delete magSfPtr;
    }
    delete magSbPtr;
}


// ************************************************************************* //
