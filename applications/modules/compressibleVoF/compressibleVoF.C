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

#include "compressibleVoF.H"
#include "localEulerDdtScheme.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(compressibleVoF, 0);
    addToRunTimeSelectionTable(solver, compressibleVoF, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::compressibleVoF::compressibleVoF(fvMesh& mesh)
:
    twoPhaseVoFSolver
    (
        mesh,
        autoPtr<twoPhaseVoFMixture>(new compressibleTwoPhaseVoFMixture(mesh))
    ),

    mixture_
    (
        refCast<compressibleTwoPhaseVoFMixture>(twoPhaseVoFSolver::mixture)
    ),

    p(mixture_.p()),

    vDot
    (
        IOobject
        (
            "vDot",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1*fvc::div(phi)
    ),

    pressureReference_
    (
        p,
        p_rgh,
        pimple.dict(),
        false
    ),

    alphaRhoPhi1
    (
        IOobject::groupName("alphaRhoPhi", alpha1.group()),
        fvc::interpolate(mixture_.thermo1().rho())*alphaPhi1
    ),

    alphaRhoPhi2
    (
        IOobject::groupName("alphaRhoPhi", alpha2.group()),
        fvc::interpolate(mixture_.thermo2().rho())*alphaPhi2
    ),

    K("K", 0.5*magSqr(U)),

    momentumTransport
    (
        rho,
        U,
        phi,
        rhoPhi,
        alphaPhi1,
        alphaPhi2,
        alphaRhoPhi1,
        alphaRhoPhi2,
        mixture_
    ),

    thermophysicalTransport(momentumTransport),

    mixture(mixture_)
{
    if (correctPhi || mesh.topoChanging())
    {
        rAU = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTime/dimDensity, 1)
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::compressibleVoF::~compressibleVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::compressibleVoF::prePredictor()
{
    twoPhaseVoFSolver::prePredictor();

    const volScalarField& rho1 = mixture_.thermo1().rho();
    const volScalarField& rho2 = mixture_.thermo2().rho();

    alphaRhoPhi1 = fvc::interpolate(rho1)*alphaPhi1;
    alphaRhoPhi2 = fvc::interpolate(rho2)*alphaPhi2;

    rhoPhi = alphaRhoPhi1 + alphaRhoPhi2;

    contErr1 =
    (
        fvc::ddt(alpha1, rho1)()() + fvc::div(alphaRhoPhi1)()()
      - (fvModels().source(alpha1, rho1)&rho1)()
    );

    contErr2 =
    (
        fvc::ddt(alpha2, rho2)()() + fvc::div(alphaRhoPhi2)()()
      - (fvModels().source(alpha2, rho2)&rho2)()
    );

    if (pimple.predictTransport())
    {
        momentumTransport.predict();
    }

    if (pimple.predictTransport())
    {
        thermophysicalTransport.predict();
    }
}


void Foam::solvers::compressibleVoF::postCorrector()
{
    if (pimple.correctTransport())
    {
        momentumTransport.correct();
        thermophysicalTransport.correct();
    }
}


// ************************************************************************* //
