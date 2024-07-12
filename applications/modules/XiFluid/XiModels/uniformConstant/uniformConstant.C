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

#include "uniformConstant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiModels
{
    defineTypeNameAndDebug(uniformConstant, 0);
    addToRunTimeSelectionTable(XiModel, uniformConstant, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::XiModels::uniformConstant::readCoeffs(const dictionary& dict)
{
    XiModel::readCoeffs(dict);

    dict.lookup("Xi") >> Xi_;
    XiModel::Xi_ == Xi_;

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiModels::uniformConstant::uniformConstant
(
    const dictionary& dict,
    const psiuMulticomponentThermo& thermo,
    const fluidThermoThermophysicalTransportModel& turbulence,
    const volScalarField& Su
)
:
    XiModel(thermo, turbulence, Su),
    Xi_(readScalar(dict.lookup("Xi")))
{
    XiModel::Xi_ == Xi_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiModels::uniformConstant::~uniformConstant()
{}


// ************************************************************************* //
