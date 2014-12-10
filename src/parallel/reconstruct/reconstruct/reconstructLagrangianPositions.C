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

#include "reconstructLagrangian.H"
#include "labelIOList.H"
#include "passiveParticleCloud.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::reconstructLagrangianPositions
(
    const polyMesh& mesh,
    const word& cloudName,
    const PtrList<fvMesh>& meshes,
    const PtrList<labelIOList>& faceProcAddressing,
    const PtrList<labelIOList>& cellProcAddressing
)
{
    passiveParticleCloud lagrangianPositions
    (
        mesh,
        cloudName,
        IDLList<passiveParticle>()
    );

    forAll(meshes, i)
    {
        const labelList& cellMap = cellProcAddressing[i];

        // faceProcAddressing not required currently.
        // const labelList& faceMap = faceProcAddressing[i];

        Cloud<passiveParticle> lpi(meshes[i], cloudName, false);

        forAllConstIter(Cloud<passiveParticle>, lpi, iter)
        {
            const passiveParticle& ppi = iter();

            // // Inverting sign if necessary and subtracting 1 from
            // // faceProcAddressing
            // label mappedTetFace = mag(faceMap[ppi.tetFace()]) - 1;

            lagrangianPositions.append
            (
                new passiveParticle
                (
                    mesh,
                    ppi.position(),
                    cellMap[ppi.cell()],
                    false
                )
            );
        }
    }

    IOPosition<Cloud<passiveParticle> >(lagrangianPositions).write();
}


// ************************************************************************* //
