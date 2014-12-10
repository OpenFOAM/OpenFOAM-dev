/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    option(name, modelType, dict, mesh),
    mode_(temperatureModeNames_.read(coeffs_.lookup("mode"))),
    Tuniform_(NULL),
    TName_("T")
{
    switch (mode_)
    {
        case tmUniform:
        {
            Tuniform_.reset
            (
                DataEntry<scalar>::New("temperature", coeffs_).ptr()
            );
            break;
        }
        case tmLookup:
        {
            TName_ = coeffs_.lookupOrDefault<word>("TName", "T");
            break;
        }
        default:
        {
            // error handling done by NamedEnum
        }
    }


    fieldNames_.setSize(1, "energy");
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::fixedTemperatureConstraint::alwaysApply() const
{
    return true;
}


void Foam::fv::fixedTemperatureConstraint::setValue
(
    fvMatrix<scalar>& eqn,
    const label
)
{
    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>("thermophysicalProperties");

    if (eqn.psi().name() == thermo.he().name())
    {
        switch (mode_)
        {
            case tmUniform:
            {
                const scalar t = mesh_.time().value();
                scalarField Tuni(cells_.size(), Tuniform_->value(t));
                eqn.setValues(cells_, thermo.he(thermo.p(), Tuni, cells_));

                break;
            }
            case tmLookup:
            {
                const volScalarField& T =
                    mesh().lookupObject<volScalarField>(TName_);

                scalarField Tlkp(T, cells_);
                eqn.setValues(cells_, thermo.he(thermo.p(), Tlkp, cells_));

                break;
            }
            default:
            {
                // error handling done by NamedEnum
            }
        }

    }
}


void Foam::fv::fixedTemperatureConstraint::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::fixedTemperatureConstraint::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        if (coeffs_.found(Tuniform_->name()))
        {
            Tuniform_.reset
            (
                DataEntry<scalar>::New(Tuniform_->name(), dict).ptr()
            );
        }

        coeffs_.readIfPresent("TName", TName_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
