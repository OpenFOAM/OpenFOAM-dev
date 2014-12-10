/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "outletMappedUniformInletFvPatchField.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::outletMappedUniformInletFvPatchField<Type>::
outletMappedUniformInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    outletPatchName_(),
    phiName_("phi")
{}


template<class Type>
Foam::outletMappedUniformInletFvPatchField<Type>::
outletMappedUniformInletFvPatchField
(
    const outletMappedUniformInletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_)
{}


template<class Type>
Foam::outletMappedUniformInletFvPatchField<Type>::
outletMappedUniformInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    outletPatchName_(dict.lookup("outletPatchName")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{}


template<class Type>
Foam::outletMappedUniformInletFvPatchField<Type>::
outletMappedUniformInletFvPatchField
(
    const outletMappedUniformInletFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_)
{}



template<class Type>
Foam::outletMappedUniformInletFvPatchField<Type>::
outletMappedUniformInletFvPatchField
(
    const outletMappedUniformInletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::outletMappedUniformInletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const GeometricField<Type, fvPatchField, volMesh>& f
    (
        dynamic_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->dimensionedInternalField()
        )
    );

    const fvPatch& p = this->patch();
    label outletPatchID =
        p.patch().boundaryMesh().findPatchID(outletPatchName_);

    if (outletPatchID < 0)
    {
        FatalErrorIn
        (
            "void outletMappedUniformInletFvPatchField<Type>::updateCoeffs()"
        )   << "Unable to find outlet patch " << outletPatchName_
            << abort(FatalError);
    }

    const fvPatch& outletPatch = p.boundaryMesh()[outletPatchID];

    const fvPatchField<Type>& outletPatchField =
        f.boundaryField()[outletPatchID];

    const surfaceScalarField& phi =
        this->db().objectRegistry::template lookupObject<surfaceScalarField>
        (phiName_);

    const scalarField& outletPatchPhi = phi.boundaryField()[outletPatchID];
    scalar sumOutletPatchPhi = gSum(outletPatchPhi);

    if (sumOutletPatchPhi > SMALL)
    {
        Type averageOutletField =
            gSum(outletPatchPhi*outletPatchField)
           /sumOutletPatchPhi;

        this->operator==(averageOutletField);
    }
    else
    {
        Type averageOutletField =
            gSum(outletPatch.magSf()*outletPatchField)
           /gSum(outletPatch.magSf());

        this->operator==(averageOutletField);
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::outletMappedUniformInletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("outletPatchName")
        << outletPatchName_ << token::END_STATEMENT << nl;
    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }
    this->writeEntry("value", os);
}


// ************************************************************************* //
