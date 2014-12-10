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

#include "atmBoundaryLayerInletEpsilonFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmBoundaryLayerInletEpsilonFvPatchScalarField::
atmBoundaryLayerInletEpsilonFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    z_(vector::zero),
    kappa_(0.41),
    Uref_(0),
    Href_(0),
    z0_(0),
    zGround_(0),
    Ustar_(0)
{}


atmBoundaryLayerInletEpsilonFvPatchScalarField::
atmBoundaryLayerInletEpsilonFvPatchScalarField
(
    const atmBoundaryLayerInletEpsilonFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    z_(ptf.z_),
    kappa_(ptf.kappa_),
    Uref_(ptf.Uref_),
    Href_(ptf.Href_),
    z0_(ptf.z0_, mapper),
    zGround_(ptf.zGround_, mapper),
    Ustar_(ptf.Ustar_, mapper)
{}


atmBoundaryLayerInletEpsilonFvPatchScalarField::
atmBoundaryLayerInletEpsilonFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    z_(dict.lookup("z")),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    Uref_(readScalar(dict.lookup("Uref"))),
    Href_(readScalar(dict.lookup("Href"))),
    z0_("z0", dict, p.size()),
    zGround_("zGround", dict, p.size()),
    Ustar_(p.size())
{
    if (mag(z_) < SMALL)
    {
        FatalErrorIn
        (
            "atmBoundaryLayerInletEpsilonFvPatchScalarField"
            "("
                "const fvPatch&, "
                "const DimensionedField<scalar, volMesh>&, "
                "const dictionary&"
            ")"
        )
            << "magnitude of z vector must be greater than zero"
            << abort(FatalError);
    }

    forAll (Ustar_, i)
    {
        Ustar_[i] = kappa_*Uref_/(log((Href_  + z0_[i])/max(z0_[i] , 0.001)));
    }

    z_ /= mag(z_);

    const vectorField& c = patch().Cf();
    scalarField::operator=(pow3(Ustar_)/(kappa_*((c & z_) - zGround_ + z0_)));
}


atmBoundaryLayerInletEpsilonFvPatchScalarField::
atmBoundaryLayerInletEpsilonFvPatchScalarField
(
    const atmBoundaryLayerInletEpsilonFvPatchScalarField& blpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(blpsf, iF),
    z_(blpsf.z_),
    kappa_(blpsf.kappa_),
    Uref_(blpsf.Uref_),
    Href_(blpsf.Href_),
    z0_(blpsf.z0_),
    zGround_(blpsf.zGround_),
    Ustar_(blpsf.Ustar_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmBoundaryLayerInletEpsilonFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    z0_.autoMap(m);
    zGround_.autoMap(m);
    Ustar_.autoMap(m);
}


void atmBoundaryLayerInletEpsilonFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const atmBoundaryLayerInletEpsilonFvPatchScalarField& blptf =
        refCast<const atmBoundaryLayerInletEpsilonFvPatchScalarField>(ptf);

    z0_.rmap(blptf.z0_, addr);
    zGround_.rmap(blptf.zGround_, addr);
    Ustar_.rmap(blptf.Ustar_, addr);
}


void atmBoundaryLayerInletEpsilonFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("z")
        << z_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa")
        << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("Uref")
        << Uref_ << token::END_STATEMENT << nl;
    os.writeKeyword("Href")
        << Href_ << token::END_STATEMENT << nl;
    z0_.writeEntry("z0", os);
    zGround_.writeEntry("zGround", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    atmBoundaryLayerInletEpsilonFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
