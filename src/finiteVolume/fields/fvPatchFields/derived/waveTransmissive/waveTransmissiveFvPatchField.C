/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "waveTransmissiveFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(p, iF),
    psiName_("thermo:psi"),
    gamma_(0.0)
{}


template<class Type>
Foam::waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const waveTransmissiveFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    advectiveFvPatchField<Type>(ptf, p, iF, mapper),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_)
{}


template<class Type>
Foam::waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    advectiveFvPatchField<Type>(p, iF, dict),
    psiName_(dict.lookupOrDefault<word>("psi", "thermo:psi")),
    gamma_(readScalar(dict.lookup("gamma")))
{}


template<class Type>
Foam::waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const waveTransmissiveFvPatchField& ptpsf
)
:
    advectiveFvPatchField<Type>(ptpsf),
    psiName_(ptpsf.psiName_),
    gamma_(ptpsf.gamma_)
{}


template<class Type>
Foam::waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const waveTransmissiveFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(ptpsf, iF),
    psiName_(ptpsf.psiName_),
    gamma_(ptpsf.gamma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::waveTransmissiveFvPatchField<Type>::advectionSpeed() const
{
    // Lookup the velocity and compressibility of the patch
    const fvPatchField<scalar>& psip =
        this->patch().template
            lookupPatchField<volScalarField, scalar>(psiName_);

    const surfaceScalarField& phi =
        this->db().template lookupObject<surfaceScalarField>(this->phiName_);

    fvsPatchField<scalar> phip =
        this->patch().template
            lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchScalarField& rhop =
            this->patch().template
                lookupPatchField<volScalarField, scalar>(this->rhoName_);

        phip /= rhop;
    }

    // Calculate the speed of the field wave w
    // by summing the component of the velocity normal to the boundary
    // and the speed of sound (sqrt(gamma_/psi)).
    return phip/this->patch().magSf() + sqrt(gamma_/psip);
}


template<class Type>
void Foam::waveTransmissiveFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    this->template
        writeEntryIfDifferent<word>(os, "phi", "phi", this->phiName_);
    this->template
        writeEntryIfDifferent<word>(os, "rho", "rho", this->rhoName_);
    this->template
        writeEntryIfDifferent<word>(os, "psi", "thermo:psi", psiName_);

    os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << nl;

    if (this->lInf_ > SMALL)
    {
        os.writeKeyword("fieldInf") << this->fieldInf_
            << token::END_STATEMENT << nl;
        os.writeKeyword("lInf") << this->lInf_
            << token::END_STATEMENT << nl;
    }

    this->writeEntry("value", os);
}


// ************************************************************************* //
