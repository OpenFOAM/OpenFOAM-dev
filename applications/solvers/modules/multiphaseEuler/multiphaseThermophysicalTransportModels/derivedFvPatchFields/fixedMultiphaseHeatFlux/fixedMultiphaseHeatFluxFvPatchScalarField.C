/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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

#include "fixedMultiphaseHeatFluxFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"

#include "phaseSystem.H"
#include "compressibleMomentumTransportModels.H"
#include "phaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedMultiphaseHeatFluxFvPatchScalarField::
fixedMultiphaseHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    q_(p.size(), 0.0),
    relax_(1.0),
    Tmin_(0.0)
{}


Foam::fixedMultiphaseHeatFluxFvPatchScalarField::
fixedMultiphaseHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    q_("q", dict, p.size()),
    relax_(dict.lookupOrDefault<scalar>("relax", 1.0)),
    Tmin_(dict.lookupOrDefault<scalar>("Tmin", 273))
{}


Foam::fixedMultiphaseHeatFluxFvPatchScalarField::
fixedMultiphaseHeatFluxFvPatchScalarField
(
    const fixedMultiphaseHeatFluxFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(psf, p, iF, mapper),
    q_(mapper(psf.q_)),
    relax_(psf.relax_),
    Tmin_(psf.Tmin_)
{}


Foam::fixedMultiphaseHeatFluxFvPatchScalarField::
fixedMultiphaseHeatFluxFvPatchScalarField
(
    const fixedMultiphaseHeatFluxFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(psf, iF),
    q_(psf.q_),
    relax_(psf.relax_),
    Tmin_(psf.Tmin_)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedMultiphaseHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Lookup the fluid model
    const phaseSystem& fluid =
            db().lookupObject<phaseSystem>(phaseSystem::propertiesName);

    const scalarField& Tp = *this;

    scalarField A(Tp.size(), scalar(0));
    scalarField B(Tp.size(), scalar(0));
    scalarField Q(Tp.size(), scalar(0));

    forAll(fluid.phases(), phasei)
    {
        const phaseModel& phase = fluid.phases()[phasei];
        const fluidThermo& thermo = phase.thermo();

        const fvPatchScalarField& alpha =
            phase.boundaryField()[patch().index()];

        const fvPatchScalarField& T =
            thermo.T().boundaryField()[patch().index()];

        const scalarField kappaEff(phase.kappaEff(patch().index()));

        if (debug)
        {
            scalarField q0(T.snGrad()*alpha*kappaEff);
            Q += q0;

            Info<< patch().name() << " " << phase.name()
                << ": Heat flux " << gMin(q0) << " - " << gMax(q0) << endl;
        }

        A += T.patchInternalField()*alpha*kappaEff*patch().deltaCoeffs();
        B += alpha*kappaEff*patch().deltaCoeffs();
    }

    if (debug)
    {
        Info<< patch().name() << " " << ": overall heat flux "
            << gMin(Q) << " - " << gMax(Q) << " W/m2, power: "
            << gSum(patch().magSf()*Q) << " W" << endl;
    }

    operator==((1 - relax_)*Tp + relax_*max(Tmin_,(q_ + A)/(B)));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::fixedMultiphaseHeatFluxFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    m(q_, q_);
}


void Foam::fixedMultiphaseHeatFluxFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const fixedMultiphaseHeatFluxFvPatchScalarField& mptf =
        refCast<const fixedMultiphaseHeatFluxFvPatchScalarField>(ptf);

    q_.rmap(mptf.q_, addr);
}


void Foam::fixedMultiphaseHeatFluxFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    fixedValueFvPatchScalarField::reset(ptf);

    const fixedMultiphaseHeatFluxFvPatchScalarField& mptf =
        refCast<const fixedMultiphaseHeatFluxFvPatchScalarField>(ptf);

    q_.reset(mptf.q_);
}


void Foam::fixedMultiphaseHeatFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntry(os, "relax", relax_);
    writeEntry(os, "q", q_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedMultiphaseHeatFluxFvPatchScalarField
    );
}


// ************************************************************************* //
