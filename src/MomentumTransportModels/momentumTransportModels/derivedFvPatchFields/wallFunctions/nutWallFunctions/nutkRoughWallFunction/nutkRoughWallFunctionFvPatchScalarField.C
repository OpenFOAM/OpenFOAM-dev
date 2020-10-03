/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
#include "momentumTransportModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

scalar nutkRoughWallFunctionFvPatchScalarField::E
(
    const scalar KsPlus,
    const scalar Cs
) const
{
    // Return fn based on non-dimensional roughness height
    if (KsPlus < 2.25)
    {
        return E_;
    }
    else if (KsPlus < 90)
    {
        return
            E_
           /pow
            (
                (KsPlus - 2.25)/87.75 + Cs*KsPlus,
                sin(0.4258*(log(KsPlus) - 0.811))
            );
    }
    else
    {
        return E_/(1 + Cs*KsPlus);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> nutkRoughWallFunctionFvPatchScalarField::nut() const
{
    const label patchi = patch().index();

    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                internalField().group()
            )
        );

    const scalarField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalarField> tnutw(new scalarField(*this));
    scalarField& nutw = tnutw.ref();

    forAll(nutw, facei)
    {
        const label celli = patch().faceCells()[facei];

        const scalar uStar = Cmu25*sqrt(k[celli]);
        const scalar KsPlus = uStar*Ks_[facei]/nuw[facei];
        const scalar E = this->E(KsPlus, Cs_[facei]);
        const scalar yPlusMin = constant::mathematical::e/E;
        const scalar yPlus = max(uStar*y[facei]/nuw[facei], yPlusMin);

        // To avoid oscillations limit the change in the wall viscosity
        nutw[facei] =
            max
            (
                min
                (
                    nuw[facei]*max(yPlus*kappa_/log(E*yPlus) - 1, 0),
                    max(2*nutw[facei], nuw[facei])
                ),
                0.5*nutw[facei]
            );

        if (debug)
        {
            Info<< "yPlus = " << yPlus
                << ", KsPlus = " << KsPlus
                << ", E = " << E
                << ", yPlusMin " << yPlusMin
                << ", yPlusLam " << yPlusLam(kappa_, E)
                << ", nutw = " << nutw[facei]
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
    const nutkRoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    Ks_(mapper(ptf.Ks_)),
    Cs_(mapper(ptf.Cs_))
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
    m(Ks_, Ks_);
    m(Cs_, Cs_);
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
    writeEntry(os, "Cs", Cs_);
    writeEntry(os, "Ks", Ks_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutkRoughWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
