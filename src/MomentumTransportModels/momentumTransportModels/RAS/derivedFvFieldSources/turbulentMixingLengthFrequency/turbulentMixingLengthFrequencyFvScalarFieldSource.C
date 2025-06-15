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

#include "turbulentMixingLengthFrequencyFvScalarFieldSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCellZone.H"
#include "volFields.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentMixingLengthFrequencyFvScalarFieldSource::
turbulentMixingLengthFrequencyFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvScalarFieldSource(iF, dict),
    mixingLength_(dict.lookup<scalar>("mixingLength", dimLength)),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09))
{}


Foam::turbulentMixingLengthFrequencyFvScalarFieldSource::
turbulentMixingLengthFrequencyFvScalarFieldSource
(
    const turbulentMixingLengthFrequencyFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvScalarFieldSource(field, iF),
    mixingLength_(field.mixingLength_),
    kName_(field.kName_),
    Cmu_(field.Cmu_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentMixingLengthFrequencyFvScalarFieldSource::
~turbulentMixingLengthFrequencyFvScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::turbulentMixingLengthFrequencyFvScalarFieldSource::sourceValue
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    const DimensionedField<scalar, volMesh> ks
    (
        this->value<scalar>(kName_, model, source)
    );

    return sqrt(ks)/(pow(Cmu_, 0.25)*mixingLength_);
}


Foam::tmp<Foam::scalarField>
Foam::turbulentMixingLengthFrequencyFvScalarFieldSource::sourceValue
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    const scalarField ks(this->value<scalar>(kName_, model, source, cells));

    return sqrt(ks)/(pow(Cmu_, 0.25)*mixingLength_);
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::turbulentMixingLengthFrequencyFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    return neg0(source);
}


Foam::tmp<Foam::scalarField>
Foam::turbulentMixingLengthFrequencyFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    return neg0(source);
}


void Foam::turbulentMixingLengthFrequencyFvScalarFieldSource::write
(
    Ostream& os
) const
{
    fvScalarFieldSource::write(os);
    writeEntry(os, "mixingLength", mixingLength_);
    writeEntryIfDifferent<word>(os, "k", "k", kName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTypeFieldSource
    (
        fvScalarFieldSource,
        turbulentMixingLengthFrequencyFvScalarFieldSource
    );
}

// ************************************************************************* //
