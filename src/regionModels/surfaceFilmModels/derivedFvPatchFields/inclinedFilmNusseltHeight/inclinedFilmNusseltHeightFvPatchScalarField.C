/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "inclinedFilmNusseltHeightFvPatchScalarField.H"
#include "volFields.H"
#include "kinematicSingleLayer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inclinedFilmNusseltHeightFvPatchScalarField::
inclinedFilmNusseltHeightFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    GammaMean_(),
    a_(),
    omega_()
{}


Foam::inclinedFilmNusseltHeightFvPatchScalarField::
inclinedFilmNusseltHeightFvPatchScalarField
(
    const inclinedFilmNusseltHeightFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    GammaMean_(ptf.GammaMean_().clone().ptr()),
    a_(ptf.a_().clone().ptr()),
    omega_(ptf.omega_().clone().ptr())
{}


Foam::inclinedFilmNusseltHeightFvPatchScalarField::
inclinedFilmNusseltHeightFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    GammaMean_(DataEntry<scalar>::New("GammaMean", dict)),
    a_(DataEntry<scalar>::New("a", dict)),
    omega_(DataEntry<scalar>::New("omega", dict))
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
}


Foam::inclinedFilmNusseltHeightFvPatchScalarField::
inclinedFilmNusseltHeightFvPatchScalarField
(
    const inclinedFilmNusseltHeightFvPatchScalarField& wmfrhpsf
)
:
    fixedValueFvPatchScalarField(wmfrhpsf),
    GammaMean_(wmfrhpsf.GammaMean_().clone().ptr()),
    a_(wmfrhpsf.a_().clone().ptr()),
    omega_(wmfrhpsf.omega_().clone().ptr())
{}


Foam::inclinedFilmNusseltHeightFvPatchScalarField::
inclinedFilmNusseltHeightFvPatchScalarField
(
    const inclinedFilmNusseltHeightFvPatchScalarField& wmfrhpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wmfrhpsf, iF),
    GammaMean_(wmfrhpsf.GammaMean_().clone().ptr()),
    a_(wmfrhpsf.a_().clone().ptr()),
    omega_(wmfrhpsf.omega_().clone().ptr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inclinedFilmNusseltHeightFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchI = patch().index();

    // retrieve the film region from the database

    const regionModels::regionModel& region =
        db().time().lookupObject<regionModels::regionModel>
        (
            "surfaceFilmProperties"
        );

    const regionModels::surfaceFilmModels::kinematicSingleLayer& film =
        dynamic_cast
        <
            const regionModels::surfaceFilmModels::kinematicSingleLayer&
        >(region);

    // calculate the vector tangential to the patch

    // note: normal pointing into the domain
    const vectorField n(-patch().nf());

    // TODO: currently re-evaluating the entire gTan field to return this patch
    const scalarField gTan(film.gTan()().boundaryField()[patchI] & n);

    if (patch().size() && (max(mag(gTan)) < SMALL))
    {
        WarningIn
        (
            "void Foam::inclinedFilmNusseltHeightFvPatchScalarField::"
            "updateCoeffs()"
        )
            << "Tangential gravity component is zero.  This boundary condition "
            << "is designed to operate on patches inclined with respect to "
            << "gravity"
            << endl;
    }

    const volVectorField& nHat = film.nHat();

    const vectorField nHatp(nHat.boundaryField()[patchI].patchInternalField());

    vectorField nTan(nHatp ^ n);
    nTan /= mag(nTan) + ROOTVSMALL;

    // calculate distance in patch tangential direction

    const vectorField& Cf = patch().Cf();
    scalarField d(nTan & Cf);

    // calculate the wavy film height

    const scalar t = db().time().timeOutputValue();

    const scalar GMean = GammaMean_->value(t);
    const scalar a = a_->value(t);
    const scalar omega = omega_->value(t);

    const scalarField G(GMean + a*sin(omega*constant::mathematical::twoPi*d));

    const volScalarField& mu = film.mu();
    const scalarField mup(mu.boundaryField()[patchI].patchInternalField());

    const volScalarField& rho = film.rho();
    const scalarField rhop(rho.boundaryField()[patchI].patchInternalField());

    const scalarField Re(max(G, scalar(0.0))/mup);

    operator==
    (
        pow(3.0*sqr(mup/rhop)/(gTan + ROOTVSMALL), 0.333)*pow(Re, 0.333)
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::inclinedFilmNusseltHeightFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedValueFvPatchScalarField::write(os);
    GammaMean_->writeData(os);
    a_->writeData(os);
    omega_->writeData(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        inclinedFilmNusseltHeightFvPatchScalarField
    );
}


// ************************************************************************* //
