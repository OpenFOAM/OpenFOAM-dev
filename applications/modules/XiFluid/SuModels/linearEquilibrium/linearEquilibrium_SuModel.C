/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "linearEquilibrium_SuModel.H"
#include "laminarFlameSpeed.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace SuModels
{
    defineTypeNameAndDebug(linearEquilibrium, 0);
    addToRunTimeSelectionTable(SuModel, linearEquilibrium, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::SuModels::linearEquilibrium::readCoeffs(const dictionary& dict)
{
    return SuModel::readCoeffs(dict);

    sigmaExt_.read(dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SuModels::linearEquilibrium::linearEquilibrium
(
    const dictionary& dict,
    const psiuMulticomponentThermo& thermo,
    const fluidThermoThermophysicalTransportModel& turbulence
)
:
    unstrained(dict, thermo, turbulence),
    sigmaExt_("sigmaExt", dimless/dimTime, dict)
{
    Su_.writeOpt() = IOobject::AUTO_WRITE;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SuModels::linearEquilibrium::~linearEquilibrium()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::SuModels::linearEquilibrium::correct()
{
    const fvMesh& mesh(thermo_.mesh());

    const volVectorField& U(turbulence_.U());
    const volVectorField& n = mesh.lookupObject<volVectorField>("n");
    const volScalarField& Xi = mesh.lookupObject<volScalarField>("Xi");

    const volScalarField sigmas
    (
        ((n & n)*fvc::div(U) - (n & fvc::grad(U) & n))/Xi
      + (
            (n & n)*fvc::div(Su_*n)
          - (n & fvc::grad(Su_*n) & n)
        )*(Xi + scalar(1))/(2*Xi)
    );

    Su_ == Su0_()()*max(scalar(1) - sigmas/sigmaExt_, scalar(0.01));
}


// ************************************************************************* //
