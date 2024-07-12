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

#include "equilibrium.H"
#include "XiEqModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiModels
{
    defineTypeNameAndDebug(equilibrium, 0);
    addToRunTimeSelectionTable(XiModel, equilibrium, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::XiModels::equilibrium::readCoeffs(const dictionary& dict)
{
    return XiModel::readCoeffs(dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiModels::equilibrium::equilibrium
(
    const dictionary& dict,
    const psiuMulticomponentThermo& thermo,
    const fluidThermoThermophysicalTransportModel& turbulence,
    const volScalarField& Su
)
:
    XiModel(thermo, turbulence, Su),
    XiEqModel_(XiEqModel::New(dict, thermo, turbulence, Su))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiModels::equilibrium::~equilibrium()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiModels::equilibrium::Db() const
{
    return XiEqModel_->Db();
}


void Foam::XiModels::equilibrium::correct()
{
    Xi_ == XiEqModel_->XiEq();
}


// ************************************************************************* //
