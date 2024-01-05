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

#include "epsilonWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::scalarField>>
Foam::epsilonWallFunctionFvPatchScalarField::calculate
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

    const scalar Cmu25 = pow025(nutw.Cmu());
    const scalar Cmu75 = pow(nutw.Cmu(), 0.75);

    // Initialise the returned field pair
    Pair<tmp<scalarField>> GandEpsilon
    (
        tmp<scalarField>(new scalarField(patch().size())),
        tmp<scalarField>(new scalarField(patch().size()))
    );
    scalarField& G = GandEpsilon.first().ref();
    scalarField& epsilon  = GandEpsilon.second().ref();

    // Calculate G and epsilon
    forAll(nutw, facei)
    {
        const label celli = patch().faceCells()[facei];

        const scalar yPlus = Cmu25*y[facei]*sqrt(k[celli])/nuw[facei];

        if (yPlus > nutw.yPlusLam())
        {
            G[facei] =
                (nutw[facei] + nuw[facei])
               *magGradUw[facei]
               *Cmu25*sqrt(k[celli])
              /(nutw.kappa()*y[facei]);

            epsilon[facei] =
                Cmu75*k[celli]*sqrt(k[celli])/(nutw.kappa()*y[facei]);
        }
        else
        {
            G[facei] = 0;

            epsilon[facei] = 2*k[celli]*nuw[facei]/sqr(y[facei]);
        }
    }

    return GandEpsilon;
}


void Foam::epsilonWallFunctionFvPatchScalarField::updateCoeffsMaster()
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
    volScalarField& epsilon =
        const_cast<volScalarField&>
        (
            static_cast<const volScalarField&>
            (
                internalField()
            )
        );

    const volScalarField::Boundary& epsilonBf = epsilon.boundaryField();

    // Make all processors build the near wall distances
    mtm.yb();

    // Evaluate all the wall functions
    PtrList<scalarField> Gpfs(epsilonBf.size());
    PtrList<scalarField> epsilonPfs(epsilonBf.size());
    forAll(epsilonBf, patchi)
    {
        if (isA<epsilonWallFunctionFvPatchScalarField>(epsilonBf[patchi]))
        {
            const epsilonWallFunctionFvPatchScalarField& epsilonPf =
                refCast<const epsilonWallFunctionFvPatchScalarField>
                (epsilonBf[patchi]);

            Pair<tmp<scalarField>> GandEpsilon = epsilonPf.calculate(mtm);

            Gpfs.set(patchi, GandEpsilon.first().ptr());
            epsilonPfs.set(patchi, GandEpsilon.second().ptr());
        }
    }

    // Average the values into the wall-adjacent cells and store
    wallCellGPtr_.reset(patchFieldsToWallCellField(Gpfs).ptr());
    wallCellEpsilonPtr_.reset(patchFieldsToWallCellField(epsilonPfs).ptr());

    // Set the fractional proportion of the value in the wall cells
    UIndirectList<scalar>(G, wallCells()) =
        (1 - wallCellFraction())*scalarField(G, wallCells())
      + wallCellFraction()*wallCellGPtr_();
    UIndirectList<scalar>(epsilon.ref(), wallCells()) =
        (1 - wallCellFraction())*scalarField(epsilon.ref(), wallCells())
      + wallCellFraction()*wallCellEpsilonPtr_();
}


void Foam::epsilonWallFunctionFvPatchScalarField::manipulateMatrixMaster
(
    fvMatrix<scalar>& matrix
)
{
    if (patch().index() != masterPatchIndex())
    {
        return;
    }

    matrix.setValues(wallCells(), wallCellEpsilonPtr_(), wallCellFraction());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::epsilonWallFunctionFvPatchScalarField::
epsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    wallCellWallFunctionFvPatchScalarField(p, iF, dict, false),
    wallCellGPtr_(nullptr),
    wallCellEpsilonPtr_(nullptr)
{
    // Apply a zero-gradient condition on start-up
    operator==(patchInternalField());
}


Foam::epsilonWallFunctionFvPatchScalarField::
epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    wallCellWallFunctionFvPatchScalarField(ptf, p, iF, mapper, true),
    wallCellGPtr_(nullptr),
    wallCellEpsilonPtr_(nullptr)
{}


Foam::epsilonWallFunctionFvPatchScalarField::
epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallCellWallFunctionFvPatchScalarField(ewfpsf, iF),
    wallCellGPtr_(nullptr),
    wallCellEpsilonPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::epsilonWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    initMaster();

    updateCoeffsMaster();

    operator==(patchInternalField());

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::epsilonWallFunctionFvPatchScalarField::manipulateMatrix
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

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        epsilonWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
