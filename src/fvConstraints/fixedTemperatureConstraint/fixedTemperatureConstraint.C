/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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

#include "fixedTemperatureConstraint.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(fixedTemperatureConstraint, 0);
        addToRunTimeSelectionTable
        (
            fvConstraint,
            fixedTemperatureConstraint,
            dictionary
        );
    }

    template<>
    const char* NamedEnum<fv::fixedTemperatureConstraint::temperatureMode, 2>::
    names[] =
    {
        "uniform",
        "lookup"
    };
}

const Foam::NamedEnum<Foam::fv::fixedTemperatureConstraint::temperatureMode, 2>
    Foam::fv::fixedTemperatureConstraint::modeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::fixedTemperatureConstraint::readCoeffs()
{
    mode_ = modeNames_.read(coeffs().lookup("mode"));

    switch (mode_)
    {
        case temperatureMode::uniform:
        {
            TValue_.reset
            (
                Function1<scalar>::New("temperature", coeffs()).ptr()
            );
            break;
        }
        case temperatureMode::lookup:
        {
            TName_ = coeffs().lookupOrDefault<word>("T", "T");
            break;
        }
    }

    phaseName_ = coeffs().lookupOrDefault<word>("phase", word::null);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fixedTemperatureConstraint::fixedTemperatureConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvConstraint(name, modelType, dict, mesh),
    set_(coeffs(), mesh),
    mode_(temperatureMode::uniform),
    TValue_(nullptr),
    TName_(word::null),
    phaseName_(word::null)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::fixedTemperatureConstraint::constrainedFields() const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>
        (
            IOobject::groupName(physicalProperties::typeName, phaseName_)
        );

    return wordList(1, thermo.he().name());
}


bool Foam::fv::fixedTemperatureConstraint::constrain
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    const labelList& cells = set_.cells();

    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>
        (
            IOobject::groupName(physicalProperties::typeName, phaseName_)
        );

    switch (mode_)
    {
        case temperatureMode::uniform:
        {
            const scalar t = mesh().time().value();
            scalarField Tuni(cells.size(), TValue_->value(t));
            eqn.setValues(cells, thermo.he(Tuni, cells));
            break;
        }
        case temperatureMode::lookup:
        {
            const volScalarField& T =
                mesh().lookupObject<volScalarField>(TName_);
            scalarField Tlkp(T, cells);
            eqn.setValues(cells, thermo.he(Tlkp, cells));
            break;
        }
    }

    return cells.size();
}


void Foam::fv::fixedTemperatureConstraint::updateMesh(const mapPolyMesh& map)
{
    set_.updateMesh(map);
}


void Foam::fv::fixedTemperatureConstraint::distribute
(
    const mapDistributePolyMesh& map
)
{
    set_.distribute(map);
}


bool Foam::fv::fixedTemperatureConstraint::movePoints()
{
    set_.movePoints();
    return true;
}


bool Foam::fv::fixedTemperatureConstraint::read(const dictionary& dict)
{
    if (fvConstraint::read(dict))
    {
        set_.read(coeffs());
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
