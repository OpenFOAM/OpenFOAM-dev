/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "waveInletOutletFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "levelSet.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "waveSuperposition.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::waveInletOutletFvPatchField<Type>::waveInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    inletValueAbove_(),
    inletValueBelow_(),
    phiName_("phi")
{}


template<class Type>
Foam::waveInletOutletFvPatchField<Type>::waveInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    inletValueAbove_(Function1<Type>::New("inletValueAbove", dict)),
    inletValueBelow_(Function1<Type>::New("inletValueBelow", dict)),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=(Field<Type>("value", dict, p.size()));
    }
    else
    {
        fvPatchField<Type>::operator=(this->patchInternalField());
    }

    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = Zero;
}


template<class Type>
Foam::waveInletOutletFvPatchField<Type>::waveInletOutletFvPatchField
(
    const waveInletOutletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    inletValueAbove_(ptf.inletValueAbove_, false),
    inletValueBelow_(ptf.inletValueBelow_, false),
    phiName_(ptf.phiName_)
{}


template<class Type>
Foam::waveInletOutletFvPatchField<Type>::waveInletOutletFvPatchField
(
    const waveInletOutletFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    inletValueAbove_(ptf.inletValueAbove_, false),
    inletValueBelow_(ptf.inletValueBelow_, false),
    phiName_(ptf.phiName_)
{}


template<class Type>
Foam::waveInletOutletFvPatchField<Type>::waveInletOutletFvPatchField
(
    const waveInletOutletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    inletValueAbove_(ptf.inletValueAbove_, false),
    inletValueBelow_(ptf.inletValueBelow_, false),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::waveInletOutletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const Field<scalar>& phip =
        this->patch().template lookupPatchField<surfaceScalarField, scalar>
        (
            phiName_
        );

    const scalar t = this->db().time().timeOutputValue();
    const waveSuperposition& waves = waveSuperposition::New(this->db());

    const pointField& localPoints = this->patch().patch().localPoints();

    this->refValue() =
        levelSetAverage
        (
            this->patch(),
            waves.height(t, this->patch().Cf()),
            waves.height(t, localPoints),
            Field<Type>(this->size(), inletValueAbove_->value(t)),
            Field<Type>(localPoints.size(), inletValueAbove_->value(t)),
            Field<Type>(this->size(), inletValueBelow_->value(t)),
            Field<Type>(localPoints.size(), inletValueBelow_->value(t))
        );

    this->valueFraction() = 1 - pos0(phip);

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::waveInletOutletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, inletValueAbove_());
    writeEntry(os, inletValueBelow_());
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
}


// ************************************************************************* //
