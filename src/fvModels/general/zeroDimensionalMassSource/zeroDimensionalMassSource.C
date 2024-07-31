/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "zeroDimensionalMassSource.H"
#include "fvCellSet.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(zeroDimensionalMassSource, 0);
    addToRunTimeSelectionTable(fvModel, zeroDimensionalMassSource, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::zeroDimensionalMassSource::readCoeffs(const dictionary& dict)
{
    zeroDimensionalMassFlowRate_.reset
    (
        Function1<scalar>::New
        (
            "massFlowRate",
            mesh().time().userUnits(),
            dimMass/dimTime,
            dict
        ).ptr()
    );
}


Foam::scalar Foam::fv::zeroDimensionalMassSource::massFlowRate() const
{
    return zeroDimensionalMassFlowRate_->value(mesh().time().value());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::zeroDimensionalMassSource::zeroDimensionalMassSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    zeroDimensionalMassSourceBase(name, modelType, mesh, dict),
    zeroDimensionalMassFlowRate_()
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::zeroDimensionalMassSource::read(const dictionary& dict)
{
    if (zeroDimensionalMassSourceBase::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
