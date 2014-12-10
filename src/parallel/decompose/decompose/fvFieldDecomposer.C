/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Foam::fvFieldDecomposer::patchFieldDecomposer::patchFieldDecomposer
(
    const labelUList& addressingSlice,
    const label addressingOffset
)
:
    directAddressing_(addressingSlice)
{
    forAll(directAddressing_, i)
    {
        // Subtract one to align addressing.
        directAddressing_[i] -= addressingOffset + 1;
    }
}


Foam::fvFieldDecomposer::processorVolPatchFieldDecomposer::
processorVolPatchFieldDecomposer
(
    const fvMesh& mesh,
    const labelUList& addressingSlice
)
:
    directAddressing_(addressingSlice.size())
{
    const labelList& own = mesh.faceOwner();
    const labelList& neighb = mesh.faceNeighbour();

    forAll(directAddressing_, i)
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
                directAddressing_[i] = neighb[ai];
            }
            else
            {
                directAddressing_[i] = own[ai];
            }
        }
        else
        {
            // This is a face that used to be on a cyclic boundary
            // but has now become a parallel patch face. I cannot
            // do the interpolation properly (I would need to look
            // up the different (face) list of data), so I will
            // just grab the value from the owner cell

            directAddressing_[i] = own[ai];
        }
    }
}


Foam::fvFieldDecomposer::processorSurfacePatchFieldDecomposer::
processorSurfacePatchFieldDecomposer
(
    const labelUList& addressingSlice
)
:
    addressing_(addressingSlice.size()),
    weights_(addressingSlice.size())
{
    forAll(addressing_, i)
    {
        addressing_[i].setSize(1);
        weights_[i].setSize(1);

        addressing_[i][0] = mag(addressingSlice[i]) - 1;
        weights_[i][0] = sign(addressingSlice[i]);
    }
}


Foam::fvFieldDecomposer::fvFieldDecomposer
(
    const fvMesh& completeMesh,
    const fvMesh& procMesh,
    const labelList& faceAddressing,
    const labelList& cellAddressing,
    const labelList& boundaryAddressing
)
:
    completeMesh_(completeMesh),
    procMesh_(procMesh),
    faceAddressing_(faceAddressing),
    cellAddressing_(cellAddressing),
    boundaryAddressing_(boundaryAddressing),
    patchFieldDecomposerPtrs_
    (
        procMesh_.boundary().size(),
        static_cast<patchFieldDecomposer*>(NULL)
    ),
    processorVolPatchFieldDecomposerPtrs_
    (
        procMesh_.boundary().size(),
        static_cast<processorVolPatchFieldDecomposer*>(NULL)
    ),
    processorSurfacePatchFieldDecomposerPtrs_
    (
        procMesh_.boundary().size(),
        static_cast<processorSurfacePatchFieldDecomposer*>(NULL)
    )
{
    forAll(boundaryAddressing_, patchi)
    {
        if
        (
            boundaryAddressing_[patchi] >= 0
        && !isA<processorLduInterface>(procMesh.boundary()[patchi])
        )
        {
            patchFieldDecomposerPtrs_[patchi] = new patchFieldDecomposer
            (
                procMesh_.boundary()[patchi].patchSlice(faceAddressing_),
                completeMesh_.boundaryMesh()
                [
                    boundaryAddressing_[patchi]
                ].start()
            );
        }
        else
        {
            processorVolPatchFieldDecomposerPtrs_[patchi] =
                new processorVolPatchFieldDecomposer
                (
                    completeMesh_,
                    procMesh_.boundary()[patchi].patchSlice(faceAddressing_)
                );

            processorSurfacePatchFieldDecomposerPtrs_[patchi] =
                new processorSurfacePatchFieldDecomposer
                (
                    static_cast<const labelUList&>
                    (
                        procMesh_.boundary()[patchi].patchSlice
                        (
                            faceAddressing_
                        )
                    )
                );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvFieldDecomposer::~fvFieldDecomposer()
{
    forAll(patchFieldDecomposerPtrs_, patchi)
    {
        if (patchFieldDecomposerPtrs_[patchi])
        {
            delete patchFieldDecomposerPtrs_[patchi];
        }
    }

    forAll(processorVolPatchFieldDecomposerPtrs_, patchi)
    {
        if (processorVolPatchFieldDecomposerPtrs_[patchi])
        {
            delete processorVolPatchFieldDecomposerPtrs_[patchi];
        }
    }

    forAll(processorSurfacePatchFieldDecomposerPtrs_, patchi)
    {
        if (processorSurfacePatchFieldDecomposerPtrs_[patchi])
        {
            delete processorSurfacePatchFieldDecomposerPtrs_[patchi];
        }
    }
}

// ************************************************************************* //
