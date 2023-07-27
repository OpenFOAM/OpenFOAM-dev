/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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

#include "diffusion.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(diffusion, 0);
    addToRunTimeSelectionTable(combustionModel, diffusion, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::diffusion::diffusion
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
:
    singleStepCombustion
    (
        modelType,
        thermo,
        turb,
        combustionProperties
    ),
    C_(this->coeffs().template lookup<scalar>("C")),
    oxidantName_(this->coeffs().template lookupOrDefault<word>("oxidant", "O2"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::diffusion::~diffusion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::diffusion::correct()
{
    this->wFuel_ ==
        dimensionedScalar(dimMass/pow3(dimLength)/dimTime, 0);

    this->fresCorrect();

    const label fuelI = this->fuelIndex();

    const volScalarField& YFuel = this->thermo().Y()[fuelI];

    if (this->thermo().containsSpecie(oxidantName_))
    {
        const volScalarField& YO2 = this->thermo().Y(oxidantName_);

        this->wFuel_ ==
            C_*this->thermo().rho()*this->turbulence().nuEff()
           *mag(fvc::grad(YFuel) & fvc::grad(YO2))
           *pos0(YFuel)*pos0(YO2);
    }
}


bool Foam::combustionModels::diffusion::read()
{
    if (singleStepCombustion::read())
    {
        this->coeffs().lookup("C") >> C_ ;
        this->coeffs().readIfPresent("oxidant", oxidantName_);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
