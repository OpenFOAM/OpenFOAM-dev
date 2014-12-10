/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "noPyrolysis.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noPyrolysis, 0);
addToRunTimeSelectionTable(pyrolysisModel, noPyrolysis, mesh);
addToRunTimeSelectionTable(pyrolysisModel, noPyrolysis, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void noPyrolysis::constructThermoChemistry()
{
    solidChemistry_.reset
    (
        basicSolidChemistryModel::New(regionMesh()).ptr()
    );

    solidThermo_.reset(&solidChemistry_->solidThermo());
    radiation_.reset(radiation::radiationModel::New(solidThermo_->T()).ptr());
}


bool noPyrolysis::read()
{
    if (pyrolysisModel::read())
    {
        // no additional info to read
        return true;
    }
    else
    {
        return false;
    }
}


bool noPyrolysis::read(const dictionary& dict)
{
    if (pyrolysisModel::read(dict))
    {
        // no additional info to read
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noPyrolysis::noPyrolysis
(
    const word& modelType,
    const fvMesh& mesh,
    const word& regionType
)
:
    pyrolysisModel(mesh, regionType),
    solidChemistry_(NULL),
    solidThermo_(NULL),
    radiation_(NULL)
{
    if (active())
    {
        constructThermoChemistry();
    }
}


noPyrolysis::noPyrolysis
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& regionType
)
:
    pyrolysisModel(mesh, regionType),
    solidChemistry_(NULL),
    solidThermo_(NULL),
    radiation_(NULL)
{
    if (active())
    {
        constructThermoChemistry();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noPyrolysis::~noPyrolysis()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void noPyrolysis::preEvolveRegion()
{
    //Do nothing
}


void noPyrolysis::evolveRegion()
{
    //Do nothing
}


const volScalarField& noPyrolysis::rho() const
{
    return solidThermo_->rho();
}


const volScalarField& noPyrolysis::T() const
{
    return solidThermo_->T();
}


const tmp<volScalarField> noPyrolysis::Cp() const
{
    return solidThermo_->Cp();
}


tmp<volScalarField> noPyrolysis::kappaRad() const
{
    return radiation_->absorptionEmission().a();
}


tmp<volScalarField> noPyrolysis::kappa() const
{
     return solidThermo_->kappa();
}


const surfaceScalarField& noPyrolysis::phiGas() const
{
    FatalErrorIn("const volScalarField& noPyrolysis::phiGas() const")
        << "phiGas field not available for " << type() << abort(FatalError);
    return surfaceScalarField::null();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
