/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2026 OpenFOAM Foundation
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

#include "massSource.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(massSource, 0);
    addToRunTimeSelectionTable(fvModel, massSource, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::massSource::readCoeffs(const dictionary& dict)
{
    zone_.read(coeffs(dict));

    massFlowRate_.reset
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::massSource::massSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    massSourceBase(name, modelType, mesh, dict),
    zone_(mesh),
    massFlowRate_()
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::cellZone& Foam::fv::massSource::zone() const
{
    return zone_.zone();
}


Foam::scalar Foam::fv::massSource::V() const
{
    return zone_.V();
}


Foam::dimensionedScalar Foam::fv::massSource::S() const
{
    return
        dimensionedScalar
        (
            dimMass/dimTime,
            massFlowRate_->value(mesh().time().value())
        );
}


bool Foam::fv::massSource::movePoints()
{
    zone_.movePoints();
    return true;
}


void Foam::fv::massSource::topoChange(const polyTopoChangeMap& map)
{
    zone_.topoChange(map);
}


void Foam::fv::massSource::mapMesh(const polyMeshMap& map)
{
    zone_.mapMesh(map);
}


void Foam::fv::massSource::distribute(const polyDistributionMap& map)
{
    zone_.distribute(map);
}


bool Foam::fv::massSource::read(const dictionary& dict)
{
    if (massSourceBase::read(dict))
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
