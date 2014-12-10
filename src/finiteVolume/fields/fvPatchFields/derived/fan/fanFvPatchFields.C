/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "fanFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Tuple2.H"
#include "polynomial.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeTemplatePatchTypeField
    (
        fvPatchScalarField,
        fanFvPatchScalarField
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<>
void Foam::fanFvPatchField<Foam::scalar>::calcFanJump()
{
    if (this->cyclicPatch().owner())
    {
        const surfaceScalarField& phi =
            db().lookupObject<surfaceScalarField>("phi");

        const fvsPatchField<scalar>& phip =
            patch().patchField<surfaceScalarField, scalar>(phi);

        scalarField Un(max(phip/patch().magSf(), scalar(0)));

        if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
        {
            Un /= patch().lookupPatchField<volScalarField, scalar>("rho");
        }

        this->jump_ = max(this->jumpTable_->value(Un), scalar(0));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<>
Foam::fanFvPatchField<Foam::scalar>::fanFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    uniformJumpFvPatchField<scalar>(p, iF)
{
    if (this->cyclicPatch().owner())
    {
        if (dict.found("f"))
        {
            // Backwards compatibility
            Istream& is = dict.lookup("f");
            is.format(IOstream::ASCII);
            scalarList f(is);

            label nPows = 0;
            forAll(f, powI)
            {
                if (mag(f[powI]) > VSMALL)
                {
                    nPows++;
                }
            }
            List<Tuple2<scalar, scalar> > coeffs(nPows);
            nPows = 0;
            forAll(f, powI)
            {
                if (mag(f[powI]) > VSMALL)
                {
                    coeffs[nPows++] = Tuple2<scalar, scalar>(f[powI], powI);
                }
            }

            this->jumpTable_.reset
            (
                new polynomial("jumpTable", coeffs)
            );
        }
        else
        {
            // Generic input constructed from dictionary
            this->jumpTable_ = DataEntry<scalar>::New("jumpTable", dict);
        }
    }

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        this->evaluate(Pstream::blocking);
    }
}


// ************************************************************************* //
