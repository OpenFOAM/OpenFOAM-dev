/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2026 OpenFOAM Foundation
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

#include "freestreamVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freestreamVelocityFvPatchVectorField::freestreamVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF, dict, false),
    dimensionedFreestreamValue_
    (
        iF.name(),
        "freestreamValue",
        p,
        iF.dimensions(),
        freestreamValue(),
        dict
    )
{
    if (dict.found("value"))
    {
        fvPatchVectorField::operator=
        (
            vectorField("value", iF.dimensions(), dict, p.size())
        );
    }
    else if (p.time().completeCase())
    {
        fvPatchVectorField::operator=(freestreamValue());
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Unable to evaluate function for incomplete case "
                "and 'value' entry not provided."
            << exit(FatalIOError);
    }

    refGrad() = Zero;
    valueFraction() = 1;
}


Foam::freestreamVelocityFvPatchVectorField::freestreamVelocityFvPatchVectorField
(
    const freestreamVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    dimensionedFreestreamValue_
    (
        ptf.dimensionedFreestreamValue_,
        p,
        freestreamValue()
    )
{}


Foam::freestreamVelocityFvPatchVectorField::freestreamVelocityFvPatchVectorField
(
    const freestreamVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, fvMesh>& iF
)
:
    mixedFvPatchVectorField(ptf, iF),
    dimensionedFreestreamValue_
    (
        ptf.dimensionedFreestreamValue_,
        freestreamValue()
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freestreamVelocityFvPatchVectorField::map
(
    const fvPatchVectorField& ptf,
    const fieldMapper& mapper
)
{
    mixedFvPatchField<vector>::map(ptf, mapper);
    dimensionedFreestreamValue_.map(!mapper.direct());
}


void Foam::freestreamVelocityFvPatchVectorField::reset
(
    const fvPatchVectorField& ptf
)
{
    mixedFvPatchField<vector>::reset(ptf);
    dimensionedFreestreamValue_.reset();
}


void Foam::freestreamVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    dimensionedFreestreamValue_.update();

    const Field<vector> Up(0.5*(patchInternalField() + field()));
    const Field<scalar> magUp(mag(Up));

    const Field<vector> nf(patch().nf());

    Field<scalar>& vf = valueFraction();

    forAll(vf, i)
    {
        if (magUp[i] > vSmall)
        {
            vf[i] = 0.5 - 0.5*(Up[i] & nf[i])/magUp[i];
        }
        else
        {
            vf[i] = 0.5;
        }
    }

    mixedFvPatchField<vector>::updateCoeffs();
}


void Foam::freestreamVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, dimensionedFreestreamValue_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        freestreamVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
