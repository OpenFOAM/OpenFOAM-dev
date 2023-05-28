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

#include "dynamicPressureFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicPressureFvPatchScalarField, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::dynamicPressureFvPatchScalarField::updateCoeffs
(
    const scalarField& p0p,
    const scalarField& K0mKp
)
{
    if (updated())
    {
        return;
    }

    if (internalField().dimensions() == dimPressure)
    {
        if (psiName_ == "none")
        {
            // Variable density and low-speed compressible flow

            const fvPatchField<scalar>& rho =
                patch().lookupPatchField<volScalarField, scalar>(rhoName_);

            operator==(p0p + rho*K0mKp);
        }
        else
        {
            // High-speed compressible flow

            const fvPatchField<scalar>& psip =
                patch().lookupPatchField<volScalarField, scalar>(psiName_);

            if (gamma_ > 1)
            {
                const scalar gM1ByG = (gamma_ - 1)/gamma_;

                operator==
                (
                    p0p/pow(scalar(1) - psip*gM1ByG*K0mKp, 1/gM1ByG)
                );
            }
            else
            {
                operator==(p0p/(scalar(1) - psip*K0mKp));
            }
        }
    }
    else if (internalField().dimensions() == dimPressure/dimDensity)
    {
        // Incompressible flow

        operator==(p0p + K0mKp);
    }
    else
    {
        FatalErrorInFunction
            << " Incorrect pressure dimensions " << internalField().dimensions()
            << nl
            << "    Should be " << dimPressure
            << " for compressible/variable density flow" << nl
            << "    or " << dimPressure/dimDensity
            << " for incompressible flow," << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicPressureFvPatchScalarField::dynamicPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    psiName_(dict.lookupOrDefault<word>("psi", "none")),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 1)),
    p0_("p0", dict, p.size())
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(p0_);
    }
}


Foam::dynamicPressureFvPatchScalarField::dynamicPressureFvPatchScalarField
(
    const dynamicPressureFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(psf, p, iF, mapper),
    rhoName_(psf.rhoName_),
    psiName_(psf.psiName_),
    gamma_(psf.gamma_),
    p0_(mapper(psf.p0_))
{}


Foam::dynamicPressureFvPatchScalarField::dynamicPressureFvPatchScalarField
(
    const dynamicPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    rhoName_(tppsf.rhoName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    p0_(tppsf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dynamicPressureFvPatchScalarField::map
(
    const fvPatchScalarField& psf,
    const fvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchScalarField::map(psf, mapper);

    const dynamicPressureFvPatchScalarField& dpsf =
        refCast<const dynamicPressureFvPatchScalarField>(psf);

    mapper(p0_, dpsf.p0_);
}


void Foam::dynamicPressureFvPatchScalarField::reset
(
    const fvPatchScalarField& psf
)
{
    fixedValueFvPatchScalarField::reset(psf);

    const dynamicPressureFvPatchScalarField& dpsf =
        refCast<const dynamicPressureFvPatchScalarField>(psf);

    p0_.reset(dpsf.p0_);
}


void Foam::dynamicPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "rho", rhoName_);
    writeEntry(os, "psi", psiName_);
    writeEntry(os, "gamma", gamma_);
    writeEntry(os, "p0", p0_);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
