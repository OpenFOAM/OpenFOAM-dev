/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "instability_XiGModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiGModels
{
    defineTypeNameAndDebug(instability, 0);
    addToRunTimeSelectionTable(XiGModel, instability, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::XiGModels::instability::readCoeffs(const dictionary& dict)
{
    XiGModel::readCoeffs(dict);

    Gin_.read(dict);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiGModels::instability::instability
(
    const dictionary& dict,
    const psiuMulticomponentThermo& thermo,
    const fluidThermoThermophysicalTransportModel& thermoTransport,
    const volScalarField& Su
)
:
    XiGModel(thermo, thermoTransport, Su),
    Gin_("Gin", dimless/dimTime, dict),
    XiGModel_(XiGModel::New(dict, thermo, thermoTransport, Su))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiGModels::instability::~instability()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiGModels::instability::G() const
{
    const volScalarField turbXiG(XiGModel_->G());
    return (Gin_*Gin_/(Gin_ + turbXiG) + turbXiG);
}


// ************************************************************************* //
