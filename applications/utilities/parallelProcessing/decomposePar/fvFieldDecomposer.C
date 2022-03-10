/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "fvFieldDecomposer.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::labelList Foam::fvFieldDecomposer::patchFieldDecomposer::alignAddressing
(
    const labelUList& addressingSlice,
    const label addressingOffset
) const
{
    labelList addressing(addressingSlice.size());

    forAll(addressing, i)
    {
        // Subtract one to align addressing.
        addressing[i] = addressingSlice[i] - (addressingOffset + 1);
    }

    return addressing;
}


Foam::fvFieldDecomposer::patchFieldDecomposer::patchFieldDecomposer
(
    const labelUList& addressingSlice,
    const label addressingOffset
)
:
    labelList(alignAddressing(addressingSlice, addressingOffset)),
    directFvPatchFieldMapper(static_cast<const labelList&>(*this))
{}


Foam::labelList Foam::fvFieldDecomposer::processorVolPatchFieldDecomposer::
alignAddressing
(
    const fvMesh& mesh,
    const labelUList& addressingSlice
) const
{
    labelList addressing(addressingSlice.size());

    const labelList& own = mesh.faceOwner();
    const labelList& neighb = mesh.faceNeighbour();

    forAll(addressing, i)
    {
        // Subtract one to align addressing.
        label ai = mag(addressingSlice[i]) - 1;

        if (ai < neighb.size())
        {
            // This is a regular face. it has been an internal face
            // of the original mesh and now it has become a face
            // on the parallel boundary.
            // Give face the value of the neighbour.

            if (addressingSlice[i] >= 0)
            {
                // I have the owner so use the neighbour value
                addressing[i] = neighb[ai];
            }
            else
            {
                addressing[i] = own[ai];
            }
        }
        else
        {
            // This is a face that used to be on a cyclic boundary
            // but has now become a parallel patch face. I cannot
            // do the interpolation properly (I would need to look
            // up the different (face) list of data), so I will
            // just grab the value from the owner cell

            addressing[i] = own[ai];
        }
    }

    return addressing;
}


Foam::fvFieldDecomposer::processorVolPatchFieldDecomposer::
processorVolPatchFieldDecomposer
(
    const fvMesh& mesh,
    const labelUList& addressingSlice
)
:
    labelList(alignAddressing(mesh, addressingSlice)),
    directFvPatchFieldMapper(static_cast<const labelList&>(*this))
{}


Foam::fvFieldDecomposer::fvFieldDecomposer
(
    const fvMesh& completeMesh,
    const fvMesh& procMesh,
    const labelList& faceAddressing,
    const labelList& cellAddressing
)
:
    completeMesh_(completeMesh),
    procMesh_(procMesh),
    faceAddressing_(faceAddressing),
    cellAddressing_(cellAddressing),
    patchFieldDecomposers_(procMesh_.boundary().size()),
    processorVolPatchFieldDecomposers_(procMesh_.boundary().size())
{
    forAll(procMesh_.boundary(), procPatchi)
    {
        const fvPatch& procPatch = procMesh.boundary()[procPatchi];

        // Determine the index of the corresponding complete patch
        label completePatchi = -1;
        if (procPatchi < completeMesh_.boundary().size())
        {
            completePatchi = procPatchi;
        }
        else if (isA<processorCyclicFvPatch>(procPatch))
        {
            completePatchi =
                refCast<const processorCyclicPolyPatch>
                (procPatch.patch()).referPatchID();
        }

        // If there is a corresponding complete patch, then create a mapper
        if (completePatchi >= 0)
        {
            patchFieldDecomposers_.set
            (
                procPatchi,
                new patchFieldDecomposer
                (
                    procPatch.patchSlice(faceAddressing_),
                    completeMesh_.boundaryMesh()[completePatchi].start()
                )
            );
        }

        // If this is a processor patch then create a processor mapper
        if (procPatchi >= completeMesh_.boundary().size())
        {
            processorVolPatchFieldDecomposers_.set
            (
                procPatchi,
                new processorVolPatchFieldDecomposer
                (
                    completeMesh_,
                    procPatch.patchSlice(faceAddressing_)
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvFieldDecomposer::~fvFieldDecomposer()
{}


// ************************************************************************* //
