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

#include "isothermalFluid.H"
#include "fvCorrectPhi.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::isothermalFluid::moveMesh()
{
    if (pimple.firstIter() || pimple.moveMeshOuterCorrectors())
    {
        // Move the mesh
        mesh_.move();

        // The rhoU and rhoU0 fields can be cleared following mesh-motion
        // now the mesh has been re-stitched as necessary
        rhoU.clear();
        rhoU0.clear();

        if (mesh.changing())
        {
            if (buoyancy.valid())
            {
                buoyancy->moveMesh();
            }

            MRF.update();

            if (correctPhi || mesh.topoChanged())
            {
                // Calculate absolute flux
                // from the mapped surface velocity
                phi = mesh.Sf() & rhoUf();

                correctUphiBCs(rho, U, phi, true);

                fv::correctPhi
                (
                    phi,
                    buoyancy.valid() ? p_rgh : p,
                    thermo.psi(),
                    autoPtr<volScalarField>(),
                    divrhoU(),
                    pimple
                );

                // Make the fluxes relative to the mesh-motion
                fvc::makeRelative(phi, rho, U);
            }

            meshCourantNo();
        }
    }
}


// ************************************************************************* //
