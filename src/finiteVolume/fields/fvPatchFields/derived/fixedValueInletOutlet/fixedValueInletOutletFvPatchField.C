/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "fixedValueInletOutletFvPatchField.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fixedValueInletOutletFvPatchField<Type>::fixedValueInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    phiName_("phi")
{}


template<class Type>
Foam::fixedValueInletOutletFvPatchField<Type>::fixedValueInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, valueRequired),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{}


template<class Type>
Foam::fixedValueInletOutletFvPatchField<Type>::fixedValueInletOutletFvPatchField
(
    const fixedValueInletOutletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper, mappingRequired),
    phiName_(ptf.phiName_)
{}


template<class Type>
Foam::fixedValueInletOutletFvPatchField<Type>::fixedValueInletOutletFvPatchField
(
    const fixedValueInletOutletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedValueInletOutletFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    // Behave as a fixed value patch where there is inflow, and fixed gradient
    // patch where there is outflow
    const scalarField& phi =
        this->patch().template
        lookupPatchField<surfaceScalarField, scalar>(phiName_);
    return (1 - pos0(phi))*Zero + pos0(phi)*pTraits<Type>::one;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedValueInletOutletFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    // Behave as a fixed value patch where there is inflow, and fixed gradient
    // patch where there is outflow
    const scalarField& phi =
        this->patch().template
        lookupPatchField<surfaceScalarField, scalar>(phiName_);
    const Field<Type> pif(this->patchInternalField());
    return (1 - pos0(phi))**this + pos0(phi)*(*this - pif);
}


template<class Type>
void Foam::fixedValueInletOutletFvPatchField<Type>::write(Ostream& os) const
{
    fixedValueFvPatchField<Type>::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
}


// ************************************************************************* //
