/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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
            option,
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
    Foam::fv::fixedTemperatureConstraint::temperatureModeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fixedTemperatureConstraint::fixedTemperatureConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    mode_(temperatureModeNames_.read(coeffs_.lookup("mode"))),
    Tuniform_(nullptr),
    TName_("T"),
    phase_(coeffs_.lookupOrDefault<word>("phase", word::null))
{
    switch (mode_)
    {
        case tmUniform:
        {
            Tuniform_.reset
            (
                Function1<scalar>::New("temperature", coeffs_).ptr()
            );
            break;
        }
        case tmLookup:
        {
            TName_ = coeffs_.lookupOrDefault<word>("T", "T");
            break;
        }
        default:
        {
            // error handling done by NamedEnum
        }
    }

    // Set the field name to that of the energy field from which the temperature
    // is obtained

    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>
        (
            IOobject::groupName(basicThermo::dictName, phase_)
        );

    fieldNames_.setSize(1, thermo.he().name());

    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::fixedTemperatureConstraint::constrain
(
    fvMatrix<scalar>& eqn,
    const label
) const
{
    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>
        (
            IOobject::groupName(basicThermo::dictName, phase_)
        );

    switch (mode_)
    {
        case tmUniform:
        {
            const scalar t = mesh_.time().value();
            scalarField Tuni(cells().size(), Tuniform_->value(t));
            eqn.setValues(cells(), thermo.he(Tuni, cells()));

            break;
        }
        case tmLookup:
        {
            const volScalarField& T =
                mesh().lookupObject<volScalarField>(TName_);

            scalarField Tlkp(T, cells());
            eqn.setValues(cells(), thermo.he(Tlkp, cells()));

            break;
        }
        default:
        {
            // error handling done by NamedEnum
        }
    }
}


bool Foam::fv::fixedTemperatureConstraint::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        if (coeffs_.found(Tuniform_->name()))
        {
            Tuniform_.reset
            (
                Function1<scalar>::New(Tuniform_->name(), dict).ptr()
            );
        }

        coeffs_.readIfPresent("T", TName_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
