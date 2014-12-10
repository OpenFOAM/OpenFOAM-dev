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

#include "inclinedFilmNusseltInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "kinematicSingleLayer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inclinedFilmNusseltInletVelocityFvPatchVectorField::
inclinedFilmNusseltInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    GammaMean_(),
    a_(),
    omega_()
{}


Foam::inclinedFilmNusseltInletVelocityFvPatchVectorField::
inclinedFilmNusseltInletVelocityFvPatchVectorField
(
    const inclinedFilmNusseltInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    GammaMean_(ptf.GammaMean_().clone().ptr()),
    a_(ptf.a_().clone().ptr()),
    omega_(ptf.omega_().clone().ptr())
{}


Foam::inclinedFilmNusseltInletVelocityFvPatchVectorField::
inclinedFilmNusseltInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    GammaMean_(DataEntry<scalar>::New("GammaMean", dict)),
    a_(DataEntry<scalar>::New("a", dict)),
    omega_(DataEntry<scalar>::New("omega", dict))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::inclinedFilmNusseltInletVelocityFvPatchVectorField::
inclinedFilmNusseltInletVelocityFvPatchVectorField
(
    const inclinedFilmNusseltInletVelocityFvPatchVectorField& fmfrpvf
)
:
    fixedValueFvPatchVectorField(fmfrpvf),
    GammaMean_(fmfrpvf.GammaMean_().clone().ptr()),
    a_(fmfrpvf.a_().clone().ptr()),
    omega_(fmfrpvf.omega_().clone().ptr())
{}


Foam::inclinedFilmNusseltInletVelocityFvPatchVectorField::
inclinedFilmNusseltInletVelocityFvPatchVectorField
(
    const inclinedFilmNusseltInletVelocityFvPatchVectorField& fmfrpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fmfrpvf, iF),
    GammaMean_(fmfrpvf.GammaMean_().clone().ptr()),
    a_(fmfrpvf.a_().clone().ptr()),
    omega_(fmfrpvf.omega_().clone().ptr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inclinedFilmNusseltInletVelocityFvPatchVectorField::updateCoeffs()
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
            "void Foam::inclinedFilmNusseltInletVelocityFvPatchVectorField::"
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

    operator==(n*pow(gTan*mup/(3.0*rhop), 0.333)*pow(Re, 0.666));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::inclinedFilmNusseltInletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
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
        fvPatchVectorField,
        inclinedFilmNusseltInletVelocityFvPatchVectorField
    );
}


// ************************************************************************* //
