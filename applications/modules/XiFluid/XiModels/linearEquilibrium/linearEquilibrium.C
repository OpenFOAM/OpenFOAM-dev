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

#include "linearEquilibrium.H"
#include "XiEqModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiModels
{
    defineTypeNameAndDebug(linearEquilibrium, 0);
    addToRunTimeSelectionTable(XiModel, linearEquilibrium, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::XiModels::linearEquilibrium::readCoeffs(const dictionary& dict)
{
    XiModel::readCoeffs(dict);

    XiShapeCoeff_ = dict.lookupOrDefault<scalar>("XiShapeCoeff", 1);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiModels::linearEquilibrium::linearEquilibrium
(
    const dictionary& dict,
    const psiuMulticomponentThermo& thermo,
    const fluidThermoThermophysicalTransportModel& turbulence,
    const volScalarField& Su
)
:
    equilibrium(dict, thermo, turbulence, Su),
    XiShapeCoeff_(dict.lookupOrDefault<scalar>("XiShapeCoeff", 1))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiModels::linearEquilibrium::~linearEquilibrium()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::XiModels::linearEquilibrium::correct()
{
    Xi_ == 1 + (1 + XiShapeCoeff_*2*(0.5 - b_))*(XiEqModel_->XiEq() - 1);
}


// ************************************************************************* //
