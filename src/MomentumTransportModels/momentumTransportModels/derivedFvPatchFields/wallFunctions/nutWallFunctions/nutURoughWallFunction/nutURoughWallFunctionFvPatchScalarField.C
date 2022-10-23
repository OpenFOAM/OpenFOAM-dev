/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "nutURoughWallFunctionFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> nutURoughWallFunctionFvPatchScalarField::nut() const
{
    const label patchi = patch().index();

    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const scalarField& y = turbModel.y()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    tmp<scalarField> tyPlus = yPlus(magUp);
    scalarField& yPlus = tyPlus.ref();

    tmp<scalarField> tnutw(new scalarField(patch().size(), 0.0));
    scalarField& nutw = tnutw.ref();

    forAll(yPlus, facei)
    {
        const scalar Re = magUp[facei]*y[facei]/nuw[facei];

        if (sqr(yPlus[facei]) > Re)
        {
            nutw[facei] = nuw[facei]*(sqr(yPlus[facei])/Re - 1);
        }
        else
        {
            nutw[facei] = 0;
        }
    }

    return tnutw;
}


tmp<scalarField> nutURoughWallFunctionFvPatchScalarField::yPlus
(
    const scalarField& magUp
) const
{
    const label patchi = patch().index();

    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const scalarField& y = turbModel.y()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    tmp<scalarField> tyPlus(new scalarField(patch().size(), 0.0));
    scalarField& yPlus = tyPlus.ref();

    static const scalar c2 = 2.25/(90 - 2.25);
    static const scalar c3 = 2.0*atan(1.0)/log(90/2.25);
    static const scalar c4 = c3*log(2.25);

    // If KsPlus is based on YPlus the extra term added to the law
    // of the wall will depend on yPlus
    forAll(yPlus, facei)
    {
        if (Ks_[facei] > 0)
        {
            // Rough Walls

            const scalar c1 = 1/(90 - 2.25) + Cs_[facei];

            const scalar magUpara = magUp[facei];
            const scalar Re = magUpara*y[facei]/nuw[facei];
            const scalar kappaRe = kappa_*Re;

            scalar yp = yPlusLam_;
            const scalar ryPlusLam = 1/yp;

            int iter = 0;
            scalar yPlusLast = 0.0;

            const scalar dKsPlusdYPlus = Ks_[facei]/y[facei];

            do
            {
                yPlusLast = yp;

                // The non-dimensional roughness height
                scalar KsPlus = yp*dKsPlusdYPlus;

                // The extra term in the law-of-the-wall
                scalar yPlusGPrime = 0;
                scalar E = E_;

                if (KsPlus >= 90)
                {
                    const scalar t1 = 1 + Cs_[facei]*KsPlus;
                    E = E_/t1;
                    yPlusGPrime = Cs_[facei]*KsPlus/t1;
                }
                else if (KsPlus > 2.25)
                {
                    const scalar t1 = c1*KsPlus - c2;
                    const scalar t2 = c3*log(KsPlus) - c4;
                    const scalar sint2 = sin(t2);
                    const scalar logt1 = log(t1);
                    E = E_/pow(t1, sint2);
                    yPlusGPrime =
                        (c1*sint2*KsPlus/t1) + (c3*logt1*cos(t2));
                }

                const scalar yPlusMin = constant::mathematical::e/E;

                if (kappa_*yPlusMin > 1)
                {
                    yp = max
                    (
                        (kappaRe + yp*(1 - yPlusGPrime))
                       /(1 + log(E*yp) - yPlusGPrime),
                        sqrt(Re)
                    );
                }
                else
                {
                    if (log(E*yp) < kappa_*yp)
                    {
                        yp = max
                        (
                            (kappaRe + yp*(1 - yPlusGPrime))
                           /(1 + log(E*yp) - yPlusGPrime),
                            yPlusMin
                        );
                    }
                    else
                    {
                        yp = sqrt(Re);
                    }
                }
            } while(mag(ryPlusLam*(yp - yPlusLast)) > 0.0001 && ++iter < 20);

            yPlus[facei] = yp;
        }
        else
        {
            // Smooth Walls
            const scalar Re = magUp[facei]*y[facei]/nuw[facei];
            const scalar ryPlusLam = 1/yPlusLam_;

            int iter = 0;
            scalar yp = yPlusLam_;
            scalar yPlusLast = yp;

            do
            {
                yPlusLast = yp;
                if (yp > yPlusLam_)
                {
                    yp = (kappa_*Re + yp)/(1 + log(E_*yp));
                }
                else
                {
                    yp = sqrt(Re);
                }
            } while(mag(ryPlusLam*(yp - yPlusLast)) > 0.0001 && ++iter < 20);

            yPlus[facei] = yp;
        }
    }

    return tyPlus;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutURoughWallFunctionFvPatchScalarField::nutURoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutUWallFunctionFvPatchScalarField(p, iF),
    Ks_(p.size(), 0.0),
    Cs_(p.size(), 0.0)
{}


nutURoughWallFunctionFvPatchScalarField::nutURoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutUWallFunctionFvPatchScalarField(p, iF, dict),
    Ks_("Ks", dict, p.size()),
    Cs_("Cs", dict, p.size())
{}


nutURoughWallFunctionFvPatchScalarField::nutURoughWallFunctionFvPatchScalarField
(
    const nutURoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutUWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    Ks_(mapper(ptf.Ks_)),
    Cs_(mapper(ptf.Cs_))
{}


nutURoughWallFunctionFvPatchScalarField::nutURoughWallFunctionFvPatchScalarField
(
    const nutURoughWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutUWallFunctionFvPatchScalarField(rwfpsf, iF),
    Ks_(rwfpsf.Ks_),
    Cs_(rwfpsf.Cs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nutURoughWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    nutUWallFunctionFvPatchScalarField::autoMap(m);
    m(Ks_, Ks_);
    m(Cs_, Cs_);
}


void nutURoughWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    nutUWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const nutURoughWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const nutURoughWallFunctionFvPatchScalarField>(ptf);

    Ks_.rmap(nrwfpsf.Ks_, addr);
    Cs_.rmap(nrwfpsf.Cs_, addr);
}


void nutURoughWallFunctionFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    nutUWallFunctionFvPatchScalarField::reset(ptf);

    const nutURoughWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const nutURoughWallFunctionFvPatchScalarField>(ptf);

    Ks_.reset(nrwfpsf.Ks_);
    Cs_.reset(nrwfpsf.Cs_);
}


void nutURoughWallFunctionFvPatchScalarField::write(Ostream& os) const
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
    nutURoughWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
