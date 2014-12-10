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

#include "pressureDirectedInletOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureDirectedInletOutletVelocityFvPatchVectorField::
pressureDirectedInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    inletDir_(p.size())
{
    refValue() = *this;
    refGrad() = vector::zero;
    valueFraction() = 0.0;
}


Foam::pressureDirectedInletOutletVelocityFvPatchVectorField::
pressureDirectedInletOutletVelocityFvPatchVectorField
(
    const pressureDirectedInletOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    inletDir_(ptf.inletDir_, mapper)
{}


Foam::pressureDirectedInletOutletVelocityFvPatchVectorField::
pressureDirectedInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    inletDir_("inletDirection", dict, p.size())
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    refValue() = *this;
    refGrad() = vector::zero;
    valueFraction() = 0.0;
}


Foam::pressureDirectedInletOutletVelocityFvPatchVectorField::
pressureDirectedInletOutletVelocityFvPatchVectorField
(
    const pressureDirectedInletOutletVelocityFvPatchVectorField& pivpvf
)
:
    mixedFvPatchVectorField(pivpvf),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    inletDir_(pivpvf.inletDir_)
{}


Foam::pressureDirectedInletOutletVelocityFvPatchVectorField::
pressureDirectedInletOutletVelocityFvPatchVectorField
(
    const pressureDirectedInletOutletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    inletDir_(pivpvf.inletDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pressureDirectedInletOutletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchVectorField::autoMap(m);
    inletDir_.autoMap(m);
}


void Foam::pressureDirectedInletOutletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    mixedFvPatchVectorField::rmap(ptf, addr);

    const pressureDirectedInletOutletVelocityFvPatchVectorField& tiptf =
        refCast<const pressureDirectedInletOutletVelocityFvPatchVectorField>
        (ptf);

    inletDir_.rmap(tiptf.inletDir_, addr);
}


void Foam::pressureDirectedInletOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    tmp<vectorField> n = patch().nf();
    tmp<scalarField> ndmagS = (n & inletDir_)*patch().magSf();

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        refValue() = inletDir_*phip/ndmagS;
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        refValue() = inletDir_*phip/(rhop*ndmagS);
    }
    else
    {
        FatalErrorIn
        (
            "pressureDirectedInletOutletVelocityFvPatchVectorField::"
            "updateCoeffs()"
        )   << "dimensions of phi are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    valueFraction() = 1.0 - pos(phip);

    mixedFvPatchVectorField::updateCoeffs();
}


void Foam::pressureDirectedInletOutletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    inletDir_.writeEntry("inletDirection", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::pressureDirectedInletOutletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=
    (
        valueFraction()*(inletDir_*(inletDir_ & pvf))
      + (1 - valueFraction())*pvf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        pressureDirectedInletOutletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
