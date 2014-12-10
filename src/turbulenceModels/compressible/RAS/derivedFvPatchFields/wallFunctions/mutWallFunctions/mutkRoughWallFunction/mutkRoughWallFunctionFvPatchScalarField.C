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

#include "mutkRoughWallFunctionFvPatchScalarField.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

scalar mutkRoughWallFunctionFvPatchScalarField::fnRough
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


tmp<scalarField> mutkRoughWallFunctionFvPatchScalarField::calcMut() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const scalarField& y = turbModel.y()[patchi];
    const scalarField& rhow = turbModel.rho().boundaryField()[patchi];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const scalarField& muw = turbModel.mu().boundaryField()[patchi];

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalarField> tmutw(new scalarField(*this));
    scalarField& mutw = tmutw();

    forAll(mutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar uStar = Cmu25*sqrt(k[faceCellI]);
        scalar yPlus = uStar*y[faceI]/(muw[faceI]/rhow[faceI]);
        scalar KsPlus = uStar*Ks_[faceI]/(muw[faceI]/rhow[faceI]);

        scalar Edash = E_;
        if (KsPlus > 2.25)
        {
            Edash /= fnRough(KsPlus, Cs_[faceI]);
        }

        scalar limitingMutw = max(mutw[faceI], muw[faceI]);

        // To avoid oscillations limit the change in the wall viscosity
        // which is particularly important if it temporarily becomes zero
        mutw[faceI] =
            max
            (
                min
                (
                    muw[faceI]
                   *(yPlus*kappa_/log(max(Edash*yPlus, 1+1e-4)) - 1),
                    2*limitingMutw
                ), 0.5*limitingMutw
            );

        if (debug)
        {
            Info<< "yPlus = " << yPlus
                << ", KsPlus = " << KsPlus
                << ", Edash = " << Edash
                << ", mutw = " << mutw[faceI]
                << endl;
        }
    }

    return tmutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mutkRoughWallFunctionFvPatchScalarField::mutkRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutkWallFunctionFvPatchScalarField(p, iF),
    Ks_(p.size(), 0.0),
    Cs_(p.size(), 0.0)
{}


mutkRoughWallFunctionFvPatchScalarField::mutkRoughWallFunctionFvPatchScalarField
(
    const mutkRoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    Ks_(ptf.Ks_, mapper),
    Cs_(ptf.Cs_, mapper)
{}


mutkRoughWallFunctionFvPatchScalarField::mutkRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mutkWallFunctionFvPatchScalarField(p, iF, dict),
    Ks_("Ks", dict, p.size()),
    Cs_("Cs", dict, p.size())
{}


mutkRoughWallFunctionFvPatchScalarField::mutkRoughWallFunctionFvPatchScalarField
(
    const mutkRoughWallFunctionFvPatchScalarField& rwfpsf
)
:
    mutkWallFunctionFvPatchScalarField(rwfpsf),
    Ks_(rwfpsf.Ks_),
    Cs_(rwfpsf.Cs_)
{}


mutkRoughWallFunctionFvPatchScalarField::mutkRoughWallFunctionFvPatchScalarField
(
    const mutkRoughWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutkWallFunctionFvPatchScalarField(rwfpsf, iF),
    Ks_(rwfpsf.Ks_),
    Cs_(rwfpsf.Cs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mutkRoughWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mutkWallFunctionFvPatchScalarField::autoMap(m);
    Ks_.autoMap(m);
    Cs_.autoMap(m);
}


void mutkRoughWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mutkWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const mutkRoughWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const mutkRoughWallFunctionFvPatchScalarField>(ptf);

    Cs_.rmap(nrwfpsf.Cs_, addr);
    Ks_.rmap(nrwfpsf.Ks_, addr);
}


void mutkRoughWallFunctionFvPatchScalarField::write(Ostream& os) const
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
    mutkRoughWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
