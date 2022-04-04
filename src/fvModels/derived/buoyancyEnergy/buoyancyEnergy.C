/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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

#include "buoyancyEnergy.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(buoyancyEnergy, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        buoyancyEnergy,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::buoyancyEnergy::readCoeffs()
{
    phaseName_ = coeffs().lookupOrDefault<word>("phase", word::null);

    UName_ =
        coeffs().lookupOrDefault<word>
        (
            "U",
            IOobject::groupName("U", phaseName_)
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::buoyancyEnergy::buoyancyEnergy
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    phaseName_(word::null),
    UName_(word::null)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::buoyancyEnergy::addSupFields() const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>
        (
            IOobject::groupName(physicalProperties::typeName, phaseName_)
        );

    return wordList(1, thermo.he().name());
}


void Foam::fv::buoyancyEnergy::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    const uniformDimensionedVectorField& g =
        mesh().lookupObject<uniformDimensionedVectorField>("g");

    const volVectorField& U = mesh().lookupObject<volVectorField>(UName_);

    eqn += rho*(U&g);
}


void Foam::fv::buoyancyEnergy::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    const uniformDimensionedVectorField& g =
        mesh().lookupObject<uniformDimensionedVectorField>("g");

    const volVectorField& U = mesh().lookupObject<volVectorField>(UName_);

    eqn += alpha*rho*(U&g);
}


bool Foam::fv::buoyancyEnergy::movePoints()
{
    return true;
}


void Foam::fv::buoyancyEnergy::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::buoyancyEnergy::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::buoyancyEnergy::distribute(const polyDistributionMap&)
{}


bool Foam::fv::buoyancyEnergy::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
