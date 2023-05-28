/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "pressureDirectedInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureDirectedInletVelocityFvPatchVectorField::
pressureDirectedInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    inletDir_("inletDirection", dict, p.size())
{}


Foam::pressureDirectedInletVelocityFvPatchVectorField::
pressureDirectedInletVelocityFvPatchVectorField
(
    const pressureDirectedInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    inletDir_(mapper(ptf.inletDir_))
{}


Foam::pressureDirectedInletVelocityFvPatchVectorField::
pressureDirectedInletVelocityFvPatchVectorField
(
    const pressureDirectedInletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    inletDir_(pivpvf.inletDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pressureDirectedInletVelocityFvPatchVectorField::map
(
    const fvPatchVectorField& ptf,
    const fvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::map(ptf, mapper);

    const pressureDirectedInletVelocityFvPatchVectorField& tiptf =
        refCast<const pressureDirectedInletVelocityFvPatchVectorField>(ptf);

    mapper(inletDir_, tiptf.inletDir_);
}


void Foam::pressureDirectedInletVelocityFvPatchVectorField::reset
(
    const fvPatchVectorField& ptf
)
{
    fixedValueFvPatchVectorField::reset(ptf);

    const pressureDirectedInletVelocityFvPatchVectorField& tiptf =
        refCast<const pressureDirectedInletVelocityFvPatchVectorField>(ptf);

    inletDir_.reset(tiptf.inletDir_);
}


void Foam::pressureDirectedInletVelocityFvPatchVectorField::updateCoeffs()
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

    if (phi.dimensions() == dimFlux)
    {
        operator==(inletDir_*phip/ndmagS);
    }
    else if (phi.dimensions() == dimMassFlux)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        operator==(inletDir_*phip/(rhop*ndmagS));
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::pressureDirectedInletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntry(os, "inletDirection", inletDir_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::pressureDirectedInletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=(inletDir_*(inletDir_ & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        pressureDirectedInletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
