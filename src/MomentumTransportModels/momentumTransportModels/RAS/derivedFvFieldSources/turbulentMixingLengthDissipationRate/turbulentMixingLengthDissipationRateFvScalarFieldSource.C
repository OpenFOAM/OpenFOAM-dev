/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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

#include "turbulentMixingLengthDissipationRateFvScalarFieldSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCellSet.H"
#include "volFields.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentMixingLengthDissipationRateFvScalarFieldSource::
turbulentMixingLengthDissipationRateFvScalarFieldSource
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


Foam::turbulentMixingLengthDissipationRateFvScalarFieldSource::
turbulentMixingLengthDissipationRateFvScalarFieldSource
(
    const turbulentMixingLengthDissipationRateFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvScalarFieldSource(field, iF),
    mixingLength_(field.mixingLength_),
    kName_(field.kName_),
    Cmu_(field.Cmu_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentMixingLengthDissipationRateFvScalarFieldSource::
~turbulentMixingLengthDissipationRateFvScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::turbulentMixingLengthDissipationRateFvScalarFieldSource::sourceValue
(
    const fvSource& source
) const
{
    const scalarField ks(this->value<scalar>(kName_, source));

    const scalar Cmu75 = pow(Cmu_, 0.75);

    return Cmu75*ks*sqrt(ks)/mixingLength_;
}


Foam::tmp<Foam::scalarField>
Foam::turbulentMixingLengthDissipationRateFvScalarFieldSource::internalCoeff
(
    const fvSource& source
) const
{
    return
        neg0(source.source(internalField().name()))
       *scalarField(source.nCells(), scalar(1));
}


void Foam::turbulentMixingLengthDissipationRateFvScalarFieldSource::write
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
        turbulentMixingLengthDissipationRateFvScalarFieldSource
    );
}

// ************************************************************************* //
