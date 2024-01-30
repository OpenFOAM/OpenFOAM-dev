/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "multiphaseEuler.H"
#include "fvcDiv.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::multiphaseEuler::moveMesh()
{
    if
    (
        pimple.flow()
     && (pimple.firstIter() || pimple.moveMeshOuterCorrectors())
    )
    {
        if
        (
            (correctPhi || mesh.topoChanged())
      // && divergent()
         && !divU.valid()
        )
        {
            // Construct and register divU for mapping
            divU = new volScalarField
            (
                "divU0",
                fvc::div
                (
                    fvc::absolute(phi, movingPhases[0].U())
                )
            );
        }

        // Move the mesh
        mesh_.move();
    }
}


void Foam::solvers::multiphaseEuler::motionCorrector()
{
    if
    (
        pimple.flow()
     && (pimple.firstIter() || pimple.moveMeshOuterCorrectors())
    )
    {
        if (mesh.changing())
        {
            buoyancy.moveMesh();

            if (correctPhi || mesh.topoChanged())
            {
                fluid_.meshUpdate();

                fluid_.correctPhi
                (
                    p_rgh,
                    divU,
                    pressureReference,
                    pimple
                );
            }

            meshCourantNo();
        }

        divU.clear();
    }
}


// ************************************************************************* //
