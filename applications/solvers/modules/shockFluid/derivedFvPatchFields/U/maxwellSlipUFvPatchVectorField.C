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

#include "maxwellSlipUFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::maxwellSlipUFvPatchVectorField::maxwellSlipUFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFixedValueSlipFvPatchVectorField(p, iF),
    TName_("T"),
    rhoName_("rho"),
    psiName_("psi"),
    muName_("mu"),
    tauMCName_("tauMC"),
    accommodationCoeff_(1.0),
    Uwall_(p.size(), vector(0.0, 0.0, 0.0)),
    thermalCreep_(true),
    curvature_(true)
{}


Foam::maxwellSlipUFvPatchVectorField::maxwellSlipUFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFixedValueSlipFvPatchVectorField(p, iF),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    psiName_(dict.lookupOrDefault<word>("psi", "psi")),
    muName_(dict.lookupOrDefault<word>("mu", "mu")),
    tauMCName_(dict.lookupOrDefault<word>("tauMC", "tauMC")),
    accommodationCoeff_(dict.lookup<scalar>("accommodationCoeff")),
    Uwall_("Uwall", dict, p.size()),
    thermalCreep_(dict.lookupOrDefault("thermalCreep", true)),
    curvature_(dict.lookupOrDefault("curvature", true))
{
    if
    (
        mag(accommodationCoeff_) < small
     || mag(accommodationCoeff_) > 2.0
    )
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "unphysical accommodationCoeff_ specified"
            << "(0 < accommodationCoeff_ <= 1)" << endl
            << exit(FatalIOError);
    }

    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );

        if (dict.found("refValue") && dict.found("valueFraction"))
        {
            this->refValue() = vectorField("refValue", dict, p.size());
            this->valueFraction() =
                scalarField("valueFraction", dict, p.size());
        }
        else
        {
            this->refValue() = *this;
            this->valueFraction() = scalar(1);
        }
    }
}


Foam::maxwellSlipUFvPatchVectorField::maxwellSlipUFvPatchVectorField
(
    const maxwellSlipUFvPatchVectorField& mspvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFixedValueSlipFvPatchVectorField(mspvf, p, iF, mapper),
    TName_(mspvf.TName_),
    rhoName_(mspvf.rhoName_),
    psiName_(mspvf.psiName_),
    muName_(mspvf.muName_),
    tauMCName_(mspvf.tauMCName_),
    accommodationCoeff_(mspvf.accommodationCoeff_),
    Uwall_(mapper(mspvf.Uwall_)),
    thermalCreep_(mspvf.thermalCreep_),
    curvature_(mspvf.curvature_)
{}


Foam::maxwellSlipUFvPatchVectorField::maxwellSlipUFvPatchVectorField
(
    const maxwellSlipUFvPatchVectorField& mspvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFixedValueSlipFvPatchVectorField(mspvf, iF),
    TName_(mspvf.TName_),
    rhoName_(mspvf.rhoName_),
    psiName_(mspvf.psiName_),
    muName_(mspvf.muName_),
    tauMCName_(mspvf.tauMCName_),
    accommodationCoeff_(mspvf.accommodationCoeff_),
    Uwall_(mspvf.Uwall_),
    thermalCreep_(mspvf.thermalCreep_),
    curvature_(mspvf.curvature_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::maxwellSlipUFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFixedValueSlipFvPatchVectorField::autoMap(m);
    m(Uwall_, Uwall_);
}


void Foam::maxwellSlipUFvPatchVectorField::rmap
(
    const fvPatchVectorField& pvf,
    const labelList& addr
)
{
    mixedFixedValueSlipFvPatchVectorField::rmap(pvf, addr);

    const maxwellSlipUFvPatchVectorField& mspvf =
        refCast<const maxwellSlipUFvPatchVectorField>(pvf);

    Uwall_.rmap(mspvf.Uwall_, addr);
}


void Foam::maxwellSlipUFvPatchVectorField::reset
(
    const fvPatchVectorField& pvf
)
{
    mixedFixedValueSlipFvPatchVectorField::reset(pvf);

    const maxwellSlipUFvPatchVectorField& mspvf =
        refCast<const maxwellSlipUFvPatchVectorField>(pvf);

    Uwall_.reset(mspvf.Uwall_);
}


void Foam::maxwellSlipUFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& pmu =
        patch().lookupPatchField<volScalarField, scalar>(muName_);
    const fvPatchScalarField& prho =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);
    const fvPatchField<scalar>& ppsi =
        patch().lookupPatchField<volScalarField, scalar>(psiName_);

    Field<scalar> C1
    (
        sqrt(ppsi*constant::mathematical::piByTwo)
      * (2.0 - accommodationCoeff_)/accommodationCoeff_
    );

    Field<scalar> pnu(pmu/prho);
    valueFraction() = (1.0/(1.0 + patch().deltaCoeffs()*C1*pnu));

    refValue() = Uwall_;

    if (thermalCreep_)
    {
        const volScalarField& vsfT =
            this->db().objectRegistry::lookupObject<volScalarField>(TName_);
        label patchi = this->patch().index();
        const fvPatchScalarField& pT = vsfT.boundaryField()[patchi];
        Field<vector> gradpT(fvc::grad(vsfT)().boundaryField()[patchi]);
        vectorField n(patch().nf());

        refValue() -= 3.0*pnu/(4.0*pT)*transform(I - n*n, gradpT);
    }

    if (curvature_)
    {
        const fvPatchTensorField& ptauMC =
            patch().lookupPatchField<volTensorField, tensor>(tauMCName_);
        vectorField n(patch().nf());

        refValue() -= C1/prho*transform(I - n*n, (n & ptauMC));
    }

    mixedFixedValueSlipFvPatchVectorField::updateCoeffs();
}


void Foam::maxwellSlipUFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "T", "T", TName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "psi", "psi", psiName_);
    writeEntryIfDifferent<word>(os, "mu", "mu", muName_);
    writeEntryIfDifferent<word>(os, "tauMC", "tauMC", tauMCName_);

    writeEntry(os, "accommodationCoeff", accommodationCoeff_);
    writeEntry(os, "Uwall", Uwall_);
    writeEntry(os, "thermalCreep", thermalCreep_);
    writeEntry(os, "curvature", curvature_);

    writeEntry(os, "refValue", refValue());
    writeEntry(os, "valueFraction", valueFraction());

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        maxwellSlipUFvPatchVectorField
    );
}

// ************************************************************************* //
