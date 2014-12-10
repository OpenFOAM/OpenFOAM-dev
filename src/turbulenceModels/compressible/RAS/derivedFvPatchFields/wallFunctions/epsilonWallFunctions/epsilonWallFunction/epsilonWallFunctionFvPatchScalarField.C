/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar epsilonWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void epsilonWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn("epsilonWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


void epsilonWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
}


void epsilonWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->dimensionedInternalField());

    const volScalarField::GeometricBoundaryField& bf = epsilon.boundaryField();

    label master = -1;
    forAll(bf, patchI)
    {
        if (isA<epsilonWallFunctionFvPatchScalarField>(bf[patchI]))
        {
            epsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchI);

            if (master == -1)
            {
                master = patchI;
            }

            epf.master() = master;
        }
    }
}


void epsilonWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->dimensionedInternalField());

    const volScalarField::GeometricBoundaryField& bf = epsilon.boundaryField();

    const fvMesh& mesh = epsilon.mesh();

    if (initialised_ && !mesh.changing())
    {
        return;
    }

    volScalarField weights
    (
        IOobject
        (
            "weights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false // do not register
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    DynamicList<label> epsilonPatches(bf.size());
    forAll(bf, patchI)
    {
        if (isA<epsilonWallFunctionFvPatchScalarField>(bf[patchI]))
        {
            epsilonPatches.append(patchI);

            const labelUList& faceCells = bf[patchI].patch().faceCells();
            forAll(faceCells, i)
            {
                label cellI = faceCells[i];
                weights[cellI]++;
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(epsilonPatches, i)
    {
        label patchI = epsilonPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchI];
        cornerWeights_[patchI] = 1.0/wf.patchInternalField();
    }

    G_.setSize(dimensionedInternalField().size(), 0.0);
    epsilon_.setSize(dimensionedInternalField().size(), 0.0);

    initialised_ = true;
}


epsilonWallFunctionFvPatchScalarField&
epsilonWallFunctionFvPatchScalarField::epsilonPatch(const label patchI)
{
    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->dimensionedInternalField());

    const volScalarField::GeometricBoundaryField& bf = epsilon.boundaryField();

    const epsilonWallFunctionFvPatchScalarField& epf =
        refCast<const epsilonWallFunctionFvPatchScalarField>(bf[patchI]);

    return const_cast<epsilonWallFunctionFvPatchScalarField&>(epf);
}


void epsilonWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbulence,
    scalarField& G0,
    scalarField& epsilon0
)
{
    // accumulate all of the G and epsilon contributions
    forAll(cornerWeights_, patchI)
    {
        if (!cornerWeights_[patchI].empty())
        {
            epsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchI);

            const List<scalar>& w = cornerWeights_[patchI];

            epf.calculate(turbulence, w, epf.patch(), G0, epsilon0);
        }
    }

    // apply zero-gradient condition for epsilon
    forAll(cornerWeights_, patchI)
    {
        if (!cornerWeights_[patchI].empty())
        {
            epsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchI);

            epf == scalarField(epsilon0, epf.patch().faceCells());
        }
    }
}


void epsilonWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbulence,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G,
    scalarField& epsilon
)
{
    const label patchI = patch.index();

    const scalarField& y = turbulence.y()[patchI];

    const scalar Cmu25 = pow025(Cmu_);
    const scalar Cmu75 = pow(Cmu_, 0.75);

    const tmp<volScalarField> tk = turbulence.k();
    const volScalarField& k = tk();

    const tmp<volScalarField> tmu = turbulence.mu();
    const scalarField& muw = tmu().boundaryField()[patchI];

    const tmp<volScalarField> tmut = turbulence.mut();
    const volScalarField& mut = tmut();
    const scalarField& mutw = mut.boundaryField()[patchI];

    const fvPatchVectorField& Uw = turbulence.U().boundaryField()[patchI];

    const scalarField magGradUw(mag(Uw.snGrad()));

    // Set epsilon and G
    forAll(mutw, faceI)
    {
        label cellI = patch.faceCells()[faceI];

        scalar w = cornerWeights[faceI];

        epsilon[cellI] += w*Cmu75*pow(k[cellI], 1.5)/(kappa_*y[faceI]);

        G[cellI] +=
            w
           *(mutw[faceI] + muw[faceI])
           *magGradUw[faceI]
           *Cmu25*sqrt(k[cellI])
           /(kappa_*y[faceI]);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    G_(),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    G_(),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    G_(),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();

    // apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ewfpsf
)
:
    fixedValueFvPatchField<scalar>(ewfpsf),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    E_(ewfpsf.E_),
    G_(),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ewfpsf, iF),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    E_(ewfpsf.E_),
    G_(),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalarField& epsilonWallFunctionFvPatchScalarField::G(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }

    return epsilonPatch(master_).G();
}


scalarField& epsilonWallFunctionFvPatchScalarField::epsilon(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            epsilon_ = 0.0;
        }

        return epsilon_;
    }

    return epsilonPatch(master_).epsilon(init);
}


void epsilonWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbulence =
        db().lookupObject<turbulenceModel>(turbulenceModel::typeName);

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbulence, G(true), epsilon(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& epsilon0 = this->epsilon();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbulence.GName())
        );

    FieldType& epsilon = const_cast<FieldType&>(dimensionedInternalField());

    forAll(*this, faceI)
    {
        label cellI = patch().faceCells()[faceI];

        G[cellI] = G0[cellI];
        epsilon[cellI] = epsilon0[cellI];
    }

    fvPatchField<scalar>::updateCoeffs();
}


void epsilonWallFunctionFvPatchScalarField::updateCoeffs
(
    const scalarField& weights
)
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbulence =
        db().lookupObject<turbulenceModel>(turbulenceModel::typeName);

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbulence, G(true), epsilon(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& epsilon0 = this->epsilon();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbulence.GName())
        );

    FieldType& epsilon = const_cast<FieldType&>(dimensionedInternalField());

    scalarField& epsilonf = *this;

    // only set the values if the weights are > tolerance
    forAll(weights, faceI)
    {
        scalar w = weights[faceI];

        if (w > tolerance_)
        {
            label cellI = patch().faceCells()[faceI];

            G[cellI] = (1.0 - w)*G[cellI] + w*G0[cellI];
            epsilon[cellI] = (1.0 - w)*epsilon[cellI] + w*epsilon0[cellI];
            epsilonf[faceI] = epsilon[cellI];
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}


void epsilonWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    matrix.setValues(patch().faceCells(), patchInternalField());

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void epsilonWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& weights
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    // filter weights so that we only apply the constraint where the
    // weight > SMALL
    DynamicList<label> constraintCells(weights.size());
    DynamicList<scalar> constraintEpsilon(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& epsilon =
        dimensionedInternalField();

    label nConstrainedCells = 0;


    forAll(weights, faceI)
    {
        // only set the values if the weights are > tolerance
        if (weights[faceI] > tolerance_)
        {
            nConstrainedCells++;

            label cellI = faceCells[faceI];

            constraintCells.append(cellI);
            constraintEpsilon.append(epsilon[cellI]);
        }
    }

    if (debug)
    {
        Pout<< "Patch: " << patch().name()
            << ": number of constrained cells = " << nConstrainedCells
            << " out of " << patch().size()
            << endl;
    }

    matrix.setValues
    (
        constraintCells,
        scalarField(constraintEpsilon.xfer())
    );

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void epsilonWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchField<scalar>::write(os);
    writeLocalEntries(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    epsilonWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
