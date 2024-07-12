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

#include "algebraic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiModels
{
    defineTypeNameAndDebug(algebraic, 0);
    addToRunTimeSelectionTable(XiModel, algebraic, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::XiModels::algebraic::readCoeffs(const dictionary& dict)
{
    XiModel::readCoeffs(dict);

    XiShapeCoeff_ = dict.lookupOrDefault<scalar>("XiShapeCoeff", 1);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiModels::algebraic::algebraic
(
    const dictionary& dict,
    const psiuMulticomponentThermo& thermo,
    const fluidThermoThermophysicalTransportModel& turbulence,
    const volScalarField& Su
)
:
    XiModel(thermo, turbulence, Su),
    XiShapeCoeff_(dict.lookupOrDefault<scalar>("XiShapeCoeff", 1)),
    XiEqModel_(XiEqModel::New(dict, thermo, turbulence, Su)),
    XiGModel_(XiGModel::New(dict, thermo, turbulence, Su))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiModels::algebraic::~algebraic()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiModels::algebraic::Db() const
{
    return XiGModel_->Db();
}


void Foam::XiModels::algebraic::correct()
{
    Xi_ == 1 + (1 + (2*XiShapeCoeff_)*(0.5 - b_))*(XiEqModel_->XiEq() - 1);
}


// ************************************************************************* //
