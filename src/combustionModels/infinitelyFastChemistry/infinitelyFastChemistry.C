/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "infinitelyFastChemistry.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(infinitelyFastChemistry, 0);
    addToRunTimeSelectionTable
    (
        combustionModel,
        infinitelyFastChemistry,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::infinitelyFastChemistry::infinitelyFastChemistry
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
    C_(this->coeffs().template lookup<scalar>("C"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::infinitelyFastChemistry::~infinitelyFastChemistry()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::infinitelyFastChemistry::correct()
{
    this->wFuel_ ==
        dimensionedScalar(dimMass/pow3(dimLength)/dimTime, 0);

    this->fresCorrect();

    const label fuelI = this->fuelIndex();

    const volScalarField& YFuel = this->thermo().Y()[fuelI];

    const dimensionedScalar s = this->s();

    if (this->thermo().containsSpecie("O2"))
    {
        const volScalarField& YO2 = this->thermo().Y("O2");

        this->wFuel_ ==
            this->rho()/(this->mesh().time().deltaT()*C_)
           *min(YFuel, YO2/s.value());
    }
}


bool Foam::combustionModels::infinitelyFastChemistry::read()
{
    if (singleStepCombustion::read())
    {
        this->coeffs().lookup("C") >> C_ ;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
