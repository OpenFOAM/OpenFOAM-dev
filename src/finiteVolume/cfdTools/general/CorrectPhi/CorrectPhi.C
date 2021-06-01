/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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

#include "CorrectPhi.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvScalarMatrix.H"
#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "fvcDiv.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "adjustPhi.H"
#include "fvcMeshPhi.H"
#include "pressureReference.H"
#include "nonOrthogonalSolutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class RAUfType, class DivUType>
void Foam::CorrectPhi
(
    surfaceScalarField& phi,
    const volVectorField& U,
    const volScalarField& p,
    const RAUfType& rAUf,
    const DivUType& divU,
    const pressureReference& pressureReference,
    nonOrthogonalSolutionControl& pcorrControl
)
{
    const fvMesh& mesh = phi.mesh();
    const Time& runTime = mesh.time();

    // Initialise BCs list for pcorr to zero-gradient
    wordList pcorrTypes
    (
        p.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    // Set BCs of pcorr to fixed-value for patches at which p is fixed
    forAll(p.boundaryField(), patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            pcorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(p.dimensions(), 0),
        pcorrTypes
    );

    if (pcorr.needReference())
    {
        fvc::makeRelative(phi, U);
        adjustPhi(phi, U, pcorr);
        fvc::makeAbsolute(phi, U);
    }

    mesh.setFluxRequired(pcorr.name());

    while (pcorrControl.correctNonOrthogonal())
    {
        // Solve for pcorr such that the divergence of the corrected flux
        // matches the divU provided (from previous iteration, time-step...)
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi) - divU
        );

        pcorrEqn.setReference(pressureReference.refCell(), 0);

        pcorrEqn.solve();

        if (pcorrControl.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }
    }
}


template<class RAUfType, class DivRhoUType>
void Foam::CorrectPhi
(
    surfaceScalarField& phi,
    const volScalarField& p,
    const volScalarField& rho,
    const volScalarField& psi,
    const RAUfType& rAUf,
    const DivRhoUType& divRhoU,
    nonOrthogonalSolutionControl& pcorrControl
)
{
    const fvMesh& mesh = phi.mesh();
    const Time& runTime = mesh.time();

    // Initialise BCs list for pcorr to zero-gradient
    wordList pcorrTypes
    (
        p.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    // Set BCs of pcorr to fixed-value for patches at which p is fixed
    forAll(p.boundaryField(), patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            pcorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(p.dimensions(), 0),
        pcorrTypes
    );

    mesh.setFluxRequired(pcorr.name());

    while (pcorrControl.correctNonOrthogonal())
    {
        // Solve for pcorr such that the divergence of the corrected flux
        // matches the divRhoU provided (from previous iteration, time-step...)
        fvScalarMatrix pcorrEqn
        (
            fvm::ddt(psi, pcorr)
          + fvc::div(phi)
          - fvm::laplacian(rAUf, pcorr)
         ==
            divRhoU
        );

        pcorrEqn.solve();

        if (pcorrControl.finalNonOrthogonalIter())
        {
            phi += pcorrEqn.flux();
        }
    }
}


// ************************************************************************* //
