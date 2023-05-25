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

#include "incompressibleFluid.H"
#include "fvCorrectPhi.H"
#include "fvcMeshPhi.H"
#include "geometricZeroField.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleFluid::moveMesh()
{
    if (pimple.firstIter() || pimple.moveMeshOuterCorrectors())
    {
        // Move the mesh
        mesh_.move();

        if (mesh.changing())
        {
            MRF.update();

            if (correctPhi || mesh.topoChanged())
            {
                // Calculate absolute flux
                // from the mapped surface velocity
                phi_ = mesh.Sf() & Uf();

                correctUphiBCs(U_, phi_, true);

                fv::correctPhi
                (
                    phi_,
                    U,
                    p,
                    autoPtr<volScalarField>(),
                    autoPtr<volScalarField>(),
                    pressureReference,
                    pimple
                );

                // Make the flux relative to the mesh motion
                fvc::makeRelative(phi_, U);
            }

            meshCourantNo();
        }
    }
}


// ************************************************************************* //
