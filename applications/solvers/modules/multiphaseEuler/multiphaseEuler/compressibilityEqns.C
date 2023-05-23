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
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSup.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::PtrList<Foam::fvScalarMatrix>
Foam::solvers::multiphaseEuler::compressibilityEqns
(
    const PtrList<volScalarField>& dmdts,
    const PtrList<volScalarField>& d2mdtdps
) const
{
    PtrList<fvScalarMatrix> pEqnComps(phases.size());

    forAll(phases_, phasei)
    {
        phaseModel& phase = phases_[phasei];
        const volScalarField& alpha = phase;
        volScalarField& rho = phase.rho();

        pEqnComps.set(phasei, new fvScalarMatrix(p_rgh, dimVolume/dimTime));
        fvScalarMatrix& pEqnComp = pEqnComps[phasei];

        // Density variation
        if (!phase.isochoric() || !phase.pure())
        {
            pEqnComp +=
                (
                    fvc::ddt(alpha, rho) + fvc::div(phase.alphaRhoPhi())
                  - fvc::Sp(fvc::ddt(alpha) + fvc::div(phase.alphaPhi()), rho)
                )/rho;
        }

        // Mesh dilatation correction
        if (mesh.moving())
        {
            pEqnComp += fvc::div(mesh.phi())*alpha;
        }

        // Compressibility
        if (!phase.incompressible())
        {
            if (pimple.transonic())
            {
                const surfaceScalarField phid
                (
                    IOobject::groupName("phid", phase.name()),
                    fvc::interpolate(phase.thermo().psi())*phase.phi()
                );

                pEqnComp +=
                    correction
                    (
                        (alpha/rho)*
                        (
                            phase.thermo().psi()*fvm::ddt(p_rgh)
                          + fvm::div(phid, p_rgh)
                          - fvm::Sp(fvc::div(phid), p_rgh)
                        )
                    );

                pEqnComps[phasei].relax();
            }
            else
            {
                pEqnComp +=
                    (alpha*phase.thermo().psi()/rho)
                   *correction(fvm::ddt(p_rgh));
            }
        }

        // Option sources
        if (fvModels().addsSupToField(rho.name()))
        {
            pEqnComp -= (fvModels().source(alpha, rho) & rho)/rho;
        }

        // Mass transfer
        if (dmdts.set(phasei))
        {
            pEqnComp -= dmdts[phasei]/rho;
        }
        if (d2mdtdps.set(phasei))
        {
            pEqnComp -= correction(fvm::Sp(d2mdtdps[phasei]/rho, p_rgh));
        }
    }

    return pEqnComps;
}


// ************************************************************************* //
