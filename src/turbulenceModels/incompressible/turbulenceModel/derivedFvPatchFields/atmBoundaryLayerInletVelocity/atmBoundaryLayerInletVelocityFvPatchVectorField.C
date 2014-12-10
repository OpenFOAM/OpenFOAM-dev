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

#include "atmBoundaryLayerInletVelocityFvPatchVectorField.H"
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

atmBoundaryLayerInletVelocityFvPatchVectorField::
atmBoundaryLayerInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Ustar_(0),
    n_(pTraits<vector>::zero),
    z_(pTraits<vector>::zero),
    z0_(0),
    kappa_(0.41),
    Uref_(0),
    Href_(0),
    zGround_(0)
{}


atmBoundaryLayerInletVelocityFvPatchVectorField::
atmBoundaryLayerInletVelocityFvPatchVectorField
(
    const atmBoundaryLayerInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Ustar_(ptf.Ustar_, mapper),
    n_(ptf.n_),
    z_(ptf.z_),
    z0_(ptf.z0_, mapper),
    kappa_(ptf.kappa_),
    Uref_(ptf.Uref_),
    Href_(ptf.Href_),
    zGround_(ptf.zGround_, mapper)
{}


atmBoundaryLayerInletVelocityFvPatchVectorField::
atmBoundaryLayerInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    Ustar_(p.size()),
    n_(dict.lookup("n")),
    z_(dict.lookup("z")),
    z0_("z0", dict, p.size()),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    Uref_(readScalar(dict.lookup("Uref"))),
    Href_(readScalar(dict.lookup("Href"))),
    zGround_("zGround", dict, p.size())
{
    if (mag(n_) < SMALL || mag(z_) < SMALL)
    {
        FatalErrorIn
        (
            "atmBoundaryLayerInletVelocityFvPatchVectorField"
            "("
                "const fvPatch&, "
                "const DimensionedField<vector, volMesh>&, "
                "onst dictionary&"
            ")"
        )
            << "magnitude of n or z must be greater than zero"
            << abort(FatalError);
    }

    n_ /= mag(n_);
    z_ /= mag(z_);

    forAll (Ustar_, i)
    {
        Ustar_[i] = kappa_*Uref_/(log((Href_  + z0_[i])/max(z0_[i] , 0.001)));
    }

    const vectorField& c = patch().Cf();
    const scalarField coord(c & z_);
    scalarField Un(coord.size());

    forAll(coord, i)
    {
        if ((coord[i] - zGround_[i]) < Href_)
        {
            Un[i] =
                (Ustar_[i]/kappa_)
              * log((coord[i] - zGround_[i] + z0_[i])/max(z0_[i], 0.001));
        }
        else
        {
            Un[i] = Uref_;
        }
    }

    vectorField::operator=(n_*Un);
}


atmBoundaryLayerInletVelocityFvPatchVectorField::
atmBoundaryLayerInletVelocityFvPatchVectorField
(
    const atmBoundaryLayerInletVelocityFvPatchVectorField& blpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(blpvf, iF),
    Ustar_(blpvf.Ustar_),
    n_(blpvf.n_),
    z_(blpvf.z_),
    z0_(blpvf.z0_),
    kappa_(blpvf.kappa_),
    Uref_(blpvf.Uref_),
    Href_(blpvf.Href_),
    zGround_(blpvf.zGround_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmBoundaryLayerInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    z0_.autoMap(m);
    zGround_.autoMap(m);
    Ustar_.autoMap(m);
}


void atmBoundaryLayerInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const atmBoundaryLayerInletVelocityFvPatchVectorField& blptf =
        refCast<const atmBoundaryLayerInletVelocityFvPatchVectorField>(ptf);

    z0_.rmap(blptf.z0_, addr);
    zGround_.rmap(blptf.zGround_, addr);
    Ustar_.rmap(blptf.Ustar_, addr);
}


void atmBoundaryLayerInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    z0_.writeEntry("z0", os) ;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("z")
        << z_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa")
        << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("Uref")
        << Uref_ << token::END_STATEMENT << nl;
    os.writeKeyword("Href")
        << Href_ << token::END_STATEMENT << nl;
    zGround_.writeEntry("zGround", os) ;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    atmBoundaryLayerInletVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
