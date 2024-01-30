/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "incompressibleDenseParticleFluid.H"
#include "fvCorrectPhi.H"
#include "fvcMeshPhi.H"
#include "geometricZeroField.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleDenseParticleFluid::moveMesh()
{
    if (pimple.firstIter() || pimple.moveMeshOuterCorrectors())
    {
        // Move the mesh
        mesh_.move();
    }
}


void Foam::solvers::incompressibleDenseParticleFluid::motionCorrector()
{
    if (pimple.firstIter() || pimple.moveMeshOuterCorrectors())
    {
        if (mesh.changing())
        {
            if (correctPhi || mesh.topoChanged())
            {
                // Calculate absolute flux
                // from the mapped surface velocity
                phic_ = mesh.Sf() & Ucf();

                correctUphiBCs(Uc_, phic_, true);

                fv::correctPhi
                (
                    phic_,
                    Uc,
                    p,
                    autoPtr<volScalarField>(),
                    autoPtr<volScalarField>(),
                    pressureReference,
                    pimple
                );

                // Make the flux relative to the mesh motion
                fvc::makeRelative(phic_, Uc);
            }

            meshCourantNo();
        }
    }
}


// ************************************************************************* //
