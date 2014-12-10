/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "nutkRoughWallFunctionFvPatchScalarField.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

scalar nutkRoughWallFunctionFvPatchScalarField::fnRough
(
    const scalar KsPlus,
    const scalar Cs
) const
{
    // Return fn based on non-dimensional roughness height

    if (KsPlus < 90.0)
    {
        return pow
        (
            (KsPlus - 2.25)/87.75 + Cs*KsPlus,
            sin(0.4258*(log(KsPlus) - 0.811))
        );
    }
    else
    {
        return (1.0 + Cs*KsPlus);
    }
}


tmp<scalarField> nutkRoughWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const scalarField& y = turbModel.y()[patchi];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const tmp<volScalarField> tnu = turbModel.nu();
    const volScalarField& nu = tnu();
    const scalarField& nuw = nu.boundaryField()[patchi];

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalarField> tnutw(new scalarField(*this));
    scalarField& nutw = tnutw();

    forAll(nutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar uStar = Cmu25*sqrt(k[faceCellI]);
        scalar yPlus = uStar*y[faceI]/nuw[faceI];
        scalar KsPlus = uStar*Ks_[faceI]/nuw[faceI];

        scalar Edash = E_;
        if (KsPlus > 2.25)
        {
            Edash /= fnRough(KsPlus, Cs_[faceI]);
        }

        scalar limitingNutw = max(nutw[faceI], nuw[faceI]);

        // To avoid oscillations limit the change in the wall viscosity
        // which is particularly important if it temporarily becomes zero
        nutw[faceI] =
            max
            (
                min
                (
                    nuw[faceI]
                   *(yPlus*kappa_/log(max(Edash*yPlus, 1+1e-4)) - 1),
                    2*limitingNutw
                ), 0.5*limitingNutw
            );

        if (debug)
        {
            Info<< "yPlus = " << yPlus
                << ", KsPlus = " << KsPlus
                << ", Edash = " << Edash
                << ", nutw = " << nutw[faceI]
                << endl;
        }
    }

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutkRoughWallFunctionFvPatchScalarField::nutkRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(p, iF),
    Ks_(p.size(), 0.0),
    Cs_(p.size(), 0.0)
{}


nutkRoughWallFunctionFvPatchScalarField::nutkRoughWallFunctionFvPatchScalarField
(
    const nutkRoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    Ks_(ptf.Ks_, mapper),
    Cs_(ptf.Cs_, mapper)
{}


nutkRoughWallFunctionFvPatchScalarField::nutkRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutkWallFunctionFvPatchScalarField(p, iF, dict),
    Ks_("Ks", dict, p.size()),
    Cs_("Cs", dict, p.size())
{}


nutkRoughWallFunctionFvPatchScalarField::nutkRoughWallFunctionFvPatchScalarField
(
    const nutkRoughWallFunctionFvPatchScalarField& rwfpsf
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf),
    Ks_(rwfpsf.Ks_),
    Cs_(rwfpsf.Cs_)
{}


nutkRoughWallFunctionFvPatchScalarField::nutkRoughWallFunctionFvPatchScalarField
(
    const nutkRoughWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf, iF),
    Ks_(rwfpsf.Ks_),
    Cs_(rwfpsf.Cs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nutkRoughWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    nutkWallFunctionFvPatchScalarField::autoMap(m);
    Ks_.autoMap(m);
    Cs_.autoMap(m);
}


void nutkRoughWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    nutkWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const nutkRoughWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const nutkRoughWallFunctionFvPatchScalarField>(ptf);

    Ks_.rmap(nrwfpsf.Ks_, addr);
    Cs_.rmap(nrwfpsf.Cs_, addr);
}


void nutkRoughWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    Cs_.writeEntry("Cs", os);
    Ks_.writeEntry("Ks", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutkRoughWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
