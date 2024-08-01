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

#include "instability.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{
    defineTypeNameAndDebug(instability, 0);
    addToRunTimeSelectionTable(XiEqModel, instability, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::XiEqModels::instability::readCoeffs(const dictionary& dict)
{
    XiEqModel::readCoeffs(dict);

    XiEqIn_.read(dict);
    lambdaIn_.read(dict);

    return XiEqModel_->read(dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiEqModels::instability::instability
(
    const dictionary& dict,
    const psiuMulticomponentThermo& thermo,
    const fluidThermoThermophysicalTransportModel& turbulence,
    const volScalarField& Su
)
:
    XiEqModel(thermo, turbulence, Su),
    XiEqIn_("XiEqIn", dimless, dict),
    lambdaIn_("lambdaIn", dimLength, dict),
    XiEqModel_(XiEqModel::New(dict, thermo, turbulence, Su))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiEqModels::instability::~instability()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiEqModels::instability::XiEq() const
{
    const volScalarField turbXiEq(XiEqModel_->XiEq());
    return XiEqIn_/turbXiEq + turbXiEq;
}


Foam::tmp<Foam::volScalarField> Foam::XiEqModels::instability::Db() const
{
    const volScalarField& rho = turbulence_.rho();

    const objectRegistry& db = Su_.db();
    const volScalarField& Xi = db.lookupObject<volScalarField>("Xi");
    const volScalarField& mgb = db.lookupObject<volScalarField>("mgb");

    return XiEqModel_->Db()
        + rho*Su_*(Xi - 1.0)*mgb*(0.5*lambdaIn_)/(mgb + 1.0/lambdaIn_);
}


// ************************************************************************* //
