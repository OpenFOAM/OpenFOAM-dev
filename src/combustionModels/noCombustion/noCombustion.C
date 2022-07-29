/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "noCombustion.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(noCombustion, 0);
    addToRunTimeSelectionTable(combustionModel, noCombustion, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::noCombustion::noCombustion
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
:
    combustionModel(modelType, thermo, turb)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::noCombustion::~noCombustion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::noCombustion::correct()
{}


Foam::tmp<Foam::fvScalarMatrix> Foam::combustionModels::noCombustion::R
(
    volScalarField& Y
) const
{
    return tmp<fvScalarMatrix>(new fvScalarMatrix(Y, dimMass/dimTime));
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::noCombustion::Qdot() const
{
    return volScalarField::New
    (
        this->thermo().phasePropertyName(typeName + ":Qdot"),
        this->mesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
    );
}


bool Foam::combustionModels::noCombustion::read()
{
    if (combustionModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
