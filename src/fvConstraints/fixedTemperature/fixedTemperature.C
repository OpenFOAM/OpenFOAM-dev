/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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

#include "fixedTemperature.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(fixedTemperature, 0);
    addToRunTimeSelectionTable
    (
        fvConstraint,
        fixedTemperature,
        dictionary
    );
    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvConstraint,
        fixedTemperature,
        dictionary,
        fixedTemperatureConstraint,
        "fixedTemperatureConstraint"
    );
}
}


namespace Foam
{
    template<>
    const char* NamedEnum<fv::fixedTemperature::temperatureMode, 2>::names[] =
        {"uniform", "lookup"};
}

const Foam::NamedEnum<Foam::fv::fixedTemperature::temperatureMode, 2>
    Foam::fv::fixedTemperature::modeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::fixedTemperature::readCoeffs()
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

    fraction_ =
        coeffs().found("fraction")
      ? Function1<scalar>::New("fraction", coeffs())
      : autoPtr<Function1<scalar>>();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fixedTemperature::fixedTemperature
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvConstraint(name, modelType, mesh, dict),
    set_(mesh, coeffs()),
    mode_(temperatureMode::uniform),
    TValue_(nullptr),
    TName_(word::null),
    phaseName_(word::null)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::fixedTemperature::constrainedFields() const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>
        (
            IOobject::groupName(physicalProperties::typeName, phaseName_)
        );

    return wordList(1, thermo.he().name());
}


bool Foam::fv::fixedTemperature::constrain
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    const labelUList cells = set_.cells();

    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>
        (
            IOobject::groupName(physicalProperties::typeName, phaseName_)
        );

    const scalar t = mesh().time().userTimeValue();

    switch (mode_)
    {
        case temperatureMode::uniform:
        {
            const scalarField Tuni(cells.size(), TValue_->value(t));
            const scalarField heuni(thermo.he(Tuni, cells));

            if (fraction_.valid())
            {
                eqn.setValues
                (
                    cells,
                    heuni,
                    scalarList(cells.size(), fraction_->value(t))
                );
            }
            else
            {
                eqn.setValues(cells, heuni);
            }

            break;
        }
        case temperatureMode::lookup:
        {
            const volScalarField& T =
                mesh().lookupObject<volScalarField>(TName_);
            const scalarField Tlkp(T, cells);
            const scalarField helkp(thermo.he(Tlkp, cells));

            if (fraction_.valid())
            {
                eqn.setValues
                (
                    cells,
                    helkp,
                    scalarList(cells.size(), fraction_->value(t))
                );
            }
            else
            {
                eqn.setValues(cells, helkp);
            }

            break;
        }
    }

    return cells.size();
}


bool Foam::fv::fixedTemperature::movePoints()
{
    set_.movePoints();
    return true;
}


void Foam::fv::fixedTemperature::topoChange(const polyTopoChangeMap& map)
{
    set_.topoChange(map);
}


void Foam::fv::fixedTemperature::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void Foam::fv::fixedTemperature::distribute(const polyDistributionMap& map)
{
    set_.distribute(map);
}


bool Foam::fv::fixedTemperature::read(const dictionary& dict)
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
