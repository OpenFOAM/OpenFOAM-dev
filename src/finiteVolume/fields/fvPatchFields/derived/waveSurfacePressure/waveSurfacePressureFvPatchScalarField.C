/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "waveSurfacePressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum
    <
        waveSurfacePressureFvPatchScalarField::timeSchemeType,
        3
    >::names[] =
    {
        fv::EulerDdtScheme<scalar>::typeName_(),
        fv::CrankNicolsonDdtScheme<scalar>::typeName_(),
        fv::backwardDdtScheme<scalar>::typeName_()
    };
}


const Foam::NamedEnum
<
    Foam::waveSurfacePressureFvPatchScalarField::timeSchemeType,
    3
>   Foam::waveSurfacePressureFvPatchScalarField::timeSchemeTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveSurfacePressureFvPatchScalarField::
waveSurfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    phiName_("phi"),
    zetaName_("zeta"),
    rhoName_("rho")
{}


Foam::waveSurfacePressureFvPatchScalarField::
waveSurfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    zetaName_(dict.lookupOrDefault<word>("zeta", "zeta")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{
    fvPatchField<scalar>::operator=
    (
        scalarField("value", dict, p.size())
    );
}


Foam::waveSurfacePressureFvPatchScalarField::
waveSurfacePressureFvPatchScalarField
(
    const waveSurfacePressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    zetaName_(ptf.zetaName_),
    rhoName_(ptf.rhoName_)
{}


Foam::waveSurfacePressureFvPatchScalarField::
waveSurfacePressureFvPatchScalarField
(
    const waveSurfacePressureFvPatchScalarField& wspsf
)
:
    fixedValueFvPatchScalarField(wspsf),
    phiName_(wspsf.phiName_),
    zetaName_(wspsf.zetaName_),
    rhoName_(wspsf.rhoName_)
{}


Foam::waveSurfacePressureFvPatchScalarField::
waveSurfacePressureFvPatchScalarField
(
    const waveSurfacePressureFvPatchScalarField& wspsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wspsf, iF),
    phiName_(wspsf.phiName_),
    zetaName_(wspsf.zetaName_),
    rhoName_(wspsf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::waveSurfacePressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchI = patch().index();

    const scalar dt = db().time().deltaTValue();

    // retrieve non-const access to zeta field from the database
    volVectorField& zeta =
        const_cast<volVectorField&>
        (
            db().lookupObject<volVectorField>(zetaName_)
        );
    vectorField& zetap = zeta.boundaryField()[patchI];

    // lookup d/dt scheme from database for zeta
    const word ddtSchemeName(zeta.mesh().ddtScheme(zeta.name()));
    timeSchemeType timeScheme(timeSchemeTypeNames_[ddtSchemeName]);

    // retrieve the flux field from the database
    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    // cache the patch face-normal vectors
    tmp<vectorField> nf(patch().nf());

    // change in zeta due to flux
    vectorField dZetap(dt*nf()*phi.boundaryField()[patchI]/patch().magSf());

    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const scalarField& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        dZetap /= rhop;
    }

    switch (timeScheme)
    {
        case tsEuler:
        case tsCrankNicolson:
        {
            zetap = zeta.oldTime().boundaryField()[patchI] + dZetap;

            break;
        }
        case tsBackward:
        {
            scalar dt0 = db().time().deltaT0Value();

            scalar c = 1.0 + dt/(dt + dt0);
            scalar c00 = dt*dt/(dt0*(dt + dt0));
            scalar c0 = c + c00;

            zetap =
                (
                    c0*zeta.oldTime().boundaryField()[patchI]
                  - c00*zeta.oldTime().oldTime().boundaryField()[patchI]
                  + dZetap
                )/c;

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "waveSurfacePressureFvPatchScalarField<Type>::updateCoeffs()"
            )   << "    Unsupported temporal differencing scheme : "
                << ddtSchemeName << nl
                << "    on patch " << this->patch().name()
                << " of field " << this->dimensionedInternalField().name()
                << " in file " << this->dimensionedInternalField().objectPath()
                << abort(FatalError);
        }
    }


    Info<< "min/max zetap = " << gMin(zetap & nf()) << ", "
        << gMax(zetap & nf()) << endl;

    // update the surface pressure
    const uniformDimensionedVectorField& g =
        db().lookupObject<uniformDimensionedVectorField>("g");

    operator==(-g.value() & zetap);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::waveSurfacePressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "zeta", "zeta", zetaName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        waveSurfacePressureFvPatchScalarField
    );
}

// ************************************************************************* //
