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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::fvFieldDecomposer::completePatchID
(
    const label procPatchi
) const
{
    const fvPatch& procPatch = procMesh_.boundary()[procPatchi];

    if (procPatchi < completeMesh_.boundary().size())
    {
        return procPatchi;
    }
    else if (isA<processorCyclicFvPatch>(procPatch))
    {
        return refCast<const processorCyclicFvPatch>(procPatch).referPatchID();
    }
    else
    {
        return -1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvFieldDecomposer::patchFieldDecomposer::patchFieldDecomposer
(
    const labelUList& addressing
)
:
    labelList(mag(addressing) - 1),
    directFvPatchFieldMapper(static_cast<const labelList&>(*this))
{}


Foam::fvFieldDecomposer::fvFieldDecomposer
(
    const fvMesh& completeMesh,
    const fvMesh& procMesh,
    const labelList& faceAddressing,
    const labelList& cellAddressing,
    const surfaceLabelField::Boundary& faceAddressingBf
)
:
    completeMesh_(completeMesh),
    procMesh_(procMesh),
    faceAddressing_(faceAddressing),
    cellAddressing_(cellAddressing),
    faceAddressingBf_(faceAddressingBf),
    patchFieldDecomposers_(procMesh_.boundary().size())
{
    forAll(procMesh_.boundary(), procPatchi)
    {
        const label completePatchi = completePatchID(procPatchi);

        // If there is a corresponding complete patch then create a patch mapper
        if (completePatchi >= 0)
        {
            patchFieldDecomposers_.set
            (
                procPatchi,
                new patchFieldDecomposer
                (
                    faceAddressingBf[completePatchi]
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvFieldDecomposer::~fvFieldDecomposer()
{}


// ************************************************************************* //
