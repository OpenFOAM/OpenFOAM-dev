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

#include "instabilityXiEq.H"
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiEqModels::instability::instability
(
    const dictionary& XiEqProperties,
    const psiuMulticomponentThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiEqModel(XiEqProperties, thermo, turbulence, Su),
    XiEqIn(XiEqModelCoeffs_.lookup<scalar>("XiEqIn")),
    XiEqModel_(XiEqModel::New(XiEqModelCoeffs_, thermo, turbulence, Su))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiEqModels::instability::~instability()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiEqModels::instability::XiEq() const
{
    volScalarField turbXiEq(XiEqModel_->XiEq());
    return XiEqIn/turbXiEq + turbXiEq;
}


bool Foam::XiEqModels::instability::read(const dictionary& XiEqProperties)
{
    XiEqModel::read(XiEqProperties);

    XiEqModelCoeffs_.lookup("XiEqIn") >> XiEqIn;

    return XiEqModel_->read(XiEqModelCoeffs_);
}


// ************************************************************************* //
