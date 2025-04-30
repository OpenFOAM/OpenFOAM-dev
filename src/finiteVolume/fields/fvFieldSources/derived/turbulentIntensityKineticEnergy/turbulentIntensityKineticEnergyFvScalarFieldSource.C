/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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

#include "turbulentIntensityKineticEnergyFvScalarFieldSource.H"
#include "fvSource.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentIntensityKineticEnergyFvScalarFieldSource::
turbulentIntensityKineticEnergyFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvScalarFieldSource(iF, dict),
    intensity_(dict.lookup<scalar>("intensity", unitFraction)),
    UName_(dict.lookupOrDefault<word>("U", "U"))
{}


Foam::turbulentIntensityKineticEnergyFvScalarFieldSource::
turbulentIntensityKineticEnergyFvScalarFieldSource
(
    const turbulentIntensityKineticEnergyFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvScalarFieldSource(field, iF),
    intensity_(field.intensity_),
    UName_(field.UName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentIntensityKineticEnergyFvScalarFieldSource::
~turbulentIntensityKineticEnergyFvScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::turbulentIntensityKineticEnergyFvScalarFieldSource::sourceValue
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    return
        1.5
       *sqr(intensity_)
       *magSqr(this->value<vector>(UName_, model, source));
}


Foam::tmp<Foam::scalarField>
Foam::turbulentIntensityKineticEnergyFvScalarFieldSource::sourceValue
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    return
        1.5
       *sqr(intensity_)
       *magSqr(this->value<vector>(UName_, model, source, cells));
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::turbulentIntensityKineticEnergyFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    return neg0(source);
}


Foam::tmp<Foam::scalarField>
Foam::turbulentIntensityKineticEnergyFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    return neg0(source);
}


void Foam::turbulentIntensityKineticEnergyFvScalarFieldSource::write
(
    Ostream& os
) const
{
    fvScalarFieldSource::write(os);
    writeEntry(os, "intensity", intensity_);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTypeFieldSource
    (
        fvScalarFieldSource,
        turbulentIntensityKineticEnergyFvScalarFieldSource
    );
}

// ************************************************************************* //
