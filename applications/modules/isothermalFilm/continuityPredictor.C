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

#include "isothermalFilm.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvmDdt.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::isothermalFilm::continuityPredictor()
{
    // Update delta and alpha BCs for time-varying inlets etc.
    delta_.correctBoundaryConditions();
    alpha_.boundaryFieldRef() == delta.boundaryField()/VbyA.boundaryField();

    fvScalarMatrix alphaEqn
    (
        fvm::ddt(rho, alpha_) + fvc::div(alphaRhoPhi)
      ==
        fvModels().source(rho, alpha_)
    );

    alphaEqn.solve();

    fvConstraints().constrain(alpha_);

    // Remove potential unboundedness in alpha caused by div(alphaRhoPhi)
    alpha_.max(0);

    // Calculate the continuity error caused by limiting alpha
    // Reset to ~0 following the alpha corrector
    correctContinuityError();

    // Update film thickness
    correctDelta();
}


void Foam::solvers::isothermalFilm::correctContinuityError()
{
    contErr =
    (
        fvc::ddt(rho, alpha)()() + fvc::div(alphaRhoPhi)()()
      - (fvModels().source(rho, alpha) & alpha)()()
    );
}


void Foam::solvers::isothermalFilm::correctDelta()
{
    delta_ = alpha*VbyA;
    delta_.correctBoundaryConditions();
}


// ************************************************************************* //
