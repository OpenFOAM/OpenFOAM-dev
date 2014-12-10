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

#include "nutkAtmRoughWallFunctionFvPatchScalarField.H"
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

tmp<scalarField> nutkAtmRoughWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchI = patch().index();

    const turbulenceModel& turbulence =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const scalarField& y = turbulence.y()[patchI];
    const tmp<volScalarField> tk = turbulence.k();
    const volScalarField& k = tk();
    const tmp<volScalarField> tnu = turbulence.nu();
    const volScalarField& nu = tnu();
    const scalarField& nuw = nu.boundaryField()[patchI];

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalarField> tnutw(new scalarField(*this));
    scalarField& nutw = tnutw();

    forAll(nutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar uStar = Cmu25*sqrt(k[faceCellI]);
        scalar yPlus = uStar*y[faceI]/nuw[faceI];

        scalar Edash = (y[faceI] + z0_[faceI])/(z0_[faceI] + 1e-4);

        nutw[faceI] =
            nuw[faceI]*(yPlus*kappa_/log(max(Edash, 1 + 1e-4)) - 1);

        if (debug)
        {
            Info<< "yPlus = " << yPlus
                << ", Edash = " << Edash
                << ", nutw = " << nutw[faceI]
                << endl;
        }
    }

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutkAtmRoughWallFunctionFvPatchScalarField::
nutkAtmRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(p, iF),
    z0_(p.size(), 0.0)
{}


nutkAtmRoughWallFunctionFvPatchScalarField::
nutkAtmRoughWallFunctionFvPatchScalarField
(
    const nutkAtmRoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    z0_(ptf.z0_, mapper)
{}


nutkAtmRoughWallFunctionFvPatchScalarField::
nutkAtmRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutkWallFunctionFvPatchScalarField(p, iF, dict),
    z0_("z0", dict, p.size())
{}


nutkAtmRoughWallFunctionFvPatchScalarField::
nutkAtmRoughWallFunctionFvPatchScalarField
(
    const nutkAtmRoughWallFunctionFvPatchScalarField& rwfpsf
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf),
    z0_(rwfpsf.z0_)
{}


nutkAtmRoughWallFunctionFvPatchScalarField::
nutkAtmRoughWallFunctionFvPatchScalarField
(
    const nutkAtmRoughWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf, iF),
    z0_(rwfpsf.z0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nutkAtmRoughWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    nutkWallFunctionFvPatchScalarField::autoMap(m);
    z0_.autoMap(m);
}


void nutkAtmRoughWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    nutkWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const nutkAtmRoughWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const nutkAtmRoughWallFunctionFvPatchScalarField>(ptf);

    z0_.rmap(nrwfpsf.z0_, addr);
}


void nutkAtmRoughWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    z0_.writeEntry("z0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutkAtmRoughWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
