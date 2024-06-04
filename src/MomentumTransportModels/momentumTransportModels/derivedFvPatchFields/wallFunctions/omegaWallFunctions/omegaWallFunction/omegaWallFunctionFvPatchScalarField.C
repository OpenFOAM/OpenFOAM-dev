/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "omegaWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::scalarField>>
Foam::omegaWallFunctionFvPatchScalarField::calculate
(
    const momentumTransportModel& mtm
) const
{
    const label patchi = patch().index();

    const nutWallFunctionFvPatchScalarField& nutw =
        nutWallFunctionFvPatchScalarField::nutw(mtm, patchi);

    const scalarField& y = mtm.yb()[patchi];

    const tmp<scalarField> tnuw = mtm.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<volScalarField> tk = mtm.k();
    const volScalarField& k = tk();

    const fvPatchVectorField& Uw = mtm.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    const VolInternalField<scalar>& Gvf =
        db().lookupObject<VolInternalField<scalar>>(mtm.GName());

    const scalar Cmu25 = pow025(nutw.Cmu());
    const scalar Cmu5 = sqrt(nutw.Cmu());

    // Initialise the returned field pair
    Pair<tmp<scalarField>> GandOmega
    (
        tmp<scalarField>(new scalarField(patch().size())),
        tmp<scalarField>(new scalarField(patch().size()))
    );
    scalarField& G = GandOmega.first().ref();
    scalarField& omega = GandOmega.second().ref();

    // Calculate G and omega
    forAll(nutw, facei)
    {
        const label celli = patch().faceCells()[facei];

        const scalar Rey = y[facei]*sqrt(k[celli])/nuw[facei];
        const scalar yPlus = Cmu25*Rey;
        const scalar uPlus = (1/nutw.kappa())*log(nutw.E()*yPlus);

        if (blended_)
        {
            const scalar lamFrac = exp(-Rey/11);
            const scalar turbFrac = 1 - lamFrac;

            const scalar uStar = sqrt
            (
                lamFrac*nuw[facei]*magGradUw[facei] + turbFrac*Cmu5*k[celli]
            );

            const scalar omegaVis = 6*nuw[facei]/(beta1_*sqr(y[facei]));
            const scalar omegaLog = uStar/(Cmu5*nutw.kappa()*y[facei]);

            G[facei] =
                   lamFrac*Gvf[celli]
                 + turbFrac
                  *sqr(uStar*magGradUw[facei]*y[facei]/uPlus)
                  /(nuw[facei]*nutw.kappa()*yPlus);

            omega[facei] = lamFrac*omegaVis + turbFrac*omegaLog;
        }
        else
        {
            if (yPlus < nutw.yPlusLam())
            {
                const scalar omegaVis = 6*nuw[facei]/(beta1_*sqr(y[facei]));

                G[facei] = Gvf[celli];

                omega[facei] = omegaVis;
            }
            else
            {
                const scalar uStar = sqrt(Cmu5*k[celli]);
                const scalar omegaLog = uStar/(Cmu5*nutw.kappa()*y[facei]);

                G[facei] =
                    sqr(uStar*magGradUw[facei]*y[facei]/uPlus)
                   /(nuw[facei]*nutw.kappa()*yPlus);

                omega[facei] = omegaLog;
            }
        }
    }

    return GandOmega;
}


void Foam::omegaWallFunctionFvPatchScalarField::updateCoeffsMaster()
{
    if (patch().index() != masterPatchIndex())
    {
        return;
    }

    // Lookup the momentum transport model
    const momentumTransportModel& mtm =
        db().lookupType<momentumTransportModel>(internalField().group());

    // Get mutable references to the turbulence fields
    VolInternalField<scalar>& G =
        db().lookupObjectRef<VolInternalField<scalar>>(mtm.GName());
    volScalarField& omega =
        const_cast<volScalarField&>
        (
            static_cast<const volScalarField&>
            (
                internalField()
            )
        );

    const volScalarField::Boundary& omegaBf = omega.boundaryField();

    // Make all processors build the near wall distances
    mtm.yb();

    // Evaluate all the wall functions
    PtrList<scalarField> Gpfs(omegaBf.size());
    PtrList<scalarField> omegaPfs(omegaBf.size());
    forAll(omegaBf, patchi)
    {
        if (isA<omegaWallFunctionFvPatchScalarField>(omegaBf[patchi]))
        {
            const omegaWallFunctionFvPatchScalarField& omegaPf =
                refCast<const omegaWallFunctionFvPatchScalarField>
                (omegaBf[patchi]);

            Pair<tmp<scalarField>> GandOmega = omegaPf.calculate(mtm);

            Gpfs.set(patchi, GandOmega.first().ptr());
            omegaPfs.set(patchi, GandOmega.second().ptr());
        }
    }

    // Average the values into the wall-adjacent cells and store
    wallCellGPtr_.reset(patchFieldsToWallCellField(Gpfs).ptr());
    wallCellOmegaPtr_.reset(patchFieldsToWallCellField(omegaPfs).ptr());

    // Set the fractional proportion of the value in the wall cells
    UIndirectList<scalar>(G, wallCells()) =
        (1 - wallCellFraction())*scalarField(G, wallCells())
      + wallCellFraction()*wallCellGPtr_();
    UIndirectList<scalar>(omega, wallCells()) =
        (1 - wallCellFraction())*scalarField(omega, wallCells())
      + wallCellFraction()*wallCellOmegaPtr_();
}


void Foam::omegaWallFunctionFvPatchScalarField::manipulateMatrixMaster
(
    fvMatrix<scalar>& matrix
)
{
    if (patch().index() != masterPatchIndex())
    {
        return;
    }

    matrix.setValues(wallCells(), wallCellOmegaPtr_(), wallCellFraction());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    wallCellWallFunctionFvPatchScalarField(p, iF, dict),
    beta1_(dict.lookupOrDefault<scalar>("beta1", 0.075)),
    blended_(dict.lookupOrDefault<Switch>("blended", false)),
    wallCellGPtr_(nullptr),
    wallCellOmegaPtr_(nullptr)
{}


Foam::omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const omegaWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    wallCellWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    beta1_(ptf.beta1_),
    blended_(ptf.blended_),
    wallCellGPtr_(nullptr),
    wallCellOmegaPtr_(nullptr)
{}


Foam::omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const omegaWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallCellWallFunctionFvPatchScalarField(owfpsf, iF),
    beta1_(owfpsf.beta1_),
    blended_(owfpsf.blended_),
    wallCellGPtr_(nullptr),
    wallCellOmegaPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::omegaWallFunctionFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    wallCellWallFunctionFvPatchScalarField::map(ptf, mapper);
    wallCellGPtr_.clear();
    wallCellOmegaPtr_.clear();
}


void Foam::omegaWallFunctionFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    wallCellWallFunctionFvPatchScalarField::reset(ptf);
    wallCellGPtr_.clear();
    wallCellOmegaPtr_.clear();
}


void Foam::omegaWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    initMaster();

    updateCoeffsMaster();

    operator==(patchInternalField());

    fvPatchScalarField::updateCoeffs();
}


void Foam::omegaWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    if (masterPatchIndex() == -1)
    {
        FatalErrorInFunction
            << "updateCoeffs must be called before manipulateMatrix"
            << exit(FatalError);
    }

    manipulateMatrixMaster(matrix);

    fvPatchScalarField::manipulateMatrix(matrix);
}


void Foam::omegaWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeEntry(os, "beta1", beta1_);
    writeEntry(os, "blended", blended_);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        omegaWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
