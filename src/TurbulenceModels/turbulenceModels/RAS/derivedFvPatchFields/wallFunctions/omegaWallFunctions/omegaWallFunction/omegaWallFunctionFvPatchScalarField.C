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

#include "omegaWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "nutkWallFunctionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar omegaWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void omegaWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn("omegaWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


void omegaWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("beta1") << beta1_ << token::END_STATEMENT << nl;
}


void omegaWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& omega =
        static_cast<const volScalarField&>(this->dimensionedInternalField());

    const volScalarField::GeometricBoundaryField& bf = omega.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<omegaWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            omegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            opf.master() = master;
        }
    }
}


void omegaWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->dimensionedInternalField());

    const volScalarField::GeometricBoundaryField& bf = omega.boundaryField();

    const fvMesh& mesh = omega.mesh();

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

    DynamicList<label> omegaPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<omegaWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            omegaPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            forAll(faceCells, i)
            {
                label cellI = faceCells[i];
                weights[cellI]++;
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(omegaPatches, i)
    {
        label patchi = omegaPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    G_.setSize(dimensionedInternalField().size(), 0.0);
    omega_.setSize(dimensionedInternalField().size(), 0.0);

    initialised_ = true;
}


omegaWallFunctionFvPatchScalarField&
omegaWallFunctionFvPatchScalarField::omegaPatch(const label patchi)
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->dimensionedInternalField());

    const volScalarField::GeometricBoundaryField& bf = omega.boundaryField();

    const omegaWallFunctionFvPatchScalarField& opf =
        refCast<const omegaWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<omegaWallFunctionFvPatchScalarField&>(opf);
}


void omegaWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbulence,
    scalarField& G0,
    scalarField& omega0
)
{
    // accumulate all of the G and omega contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            omegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            opf.calculate(turbulence, w, opf.patch(), G0, omega0);
        }
    }

    // apply zero-gradient condition for omega
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            omegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            opf == scalarField(omega0, opf.patch().faceCells());
        }
    }
}


void omegaWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbulence,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G,
    scalarField& omega
)
{
    const label patchi = patch.index();

    const scalarField& y = turbulence.y()[patchi];

    const scalar Cmu25 = pow025(Cmu_);

    const tmp<volScalarField> tk = turbulence.k();
    const volScalarField& k = tk();

    const tmp<scalarField> tnuw = turbulence.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<scalarField> tnutw = turbulence.nut(patchi);
    const scalarField& nutw = tnutw();

    const fvPatchVectorField& Uw = turbulence.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    // Set omega and G
    forAll(nutw, faceI)
    {
        label cellI = patch.faceCells()[faceI];

        scalar w = cornerWeights[faceI];

        scalar omegaVis = 6.0*nuw[faceI]/(beta1_*sqr(y[faceI]));

        scalar omegaLog = sqrt(k[cellI])/(Cmu25*kappa_*y[faceI]);

        omega[cellI] += w*sqrt(sqr(omegaVis) + sqr(omegaLog));

        G[cellI] +=
            w
           *(nutw[faceI] + nuw[faceI])
           *magGradUw[faceI]
           *Cmu25*sqrt(k[cellI])
           /(kappa_*y[faceI]);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    beta1_(0.075),
    yPlusLam_(nutkWallFunctionFvPatchScalarField::yPlusLam(kappa_, E_)),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const omegaWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    beta1_(ptf.beta1_),
    yPlusLam_(ptf.yPlusLam_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
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
    beta1_(dict.lookupOrDefault<scalar>("beta1", 0.075)),
    yPlusLam_(nutkWallFunctionFvPatchScalarField::yPlusLam(kappa_, E_)),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();

    // apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const omegaWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedValueFvPatchField<scalar>(owfpsf),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_),
    yPlusLam_(owfpsf.yPlusLam_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const omegaWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(owfpsf, iF),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_),
    yPlusLam_(owfpsf.yPlusLam_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalarField& omegaWallFunctionFvPatchScalarField::G(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }

    return omegaPatch(master_).G();
}


scalarField& omegaWallFunctionFvPatchScalarField::omega(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            omega_ = 0.0;
        }

        return omega_;
    }

    return omegaPatch(master_).omega(init);
}


void omegaWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbulence = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            dimensionedInternalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbulence, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbulence.GName())
        );

    FieldType& omega = const_cast<FieldType&>(dimensionedInternalField());

    forAll(*this, faceI)
    {
        label cellI = patch().faceCells()[faceI];

        G[cellI] = G0[cellI];
        omega[cellI] = omega0[cellI];
    }

    fvPatchField<scalar>::updateCoeffs();
}


void omegaWallFunctionFvPatchScalarField::updateCoeffs
(
    const scalarField& weights
)
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbulence = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            dimensionedInternalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbulence, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbulence.GName())
        );

    FieldType& omega = const_cast<FieldType&>(dimensionedInternalField());

    scalarField& omegaf = *this;

    // only set the values if the weights are > tolerance
    forAll(weights, faceI)
    {
        scalar w = weights[faceI];

        if (w > tolerance_)
        {
            label cellI = patch().faceCells()[faceI];

            G[cellI] = (1.0 - w)*G[cellI] + w*G0[cellI];
            omega[cellI] = (1.0 - w)*omega[cellI] + w*omega0[cellI];
            omegaf[faceI] = omega[cellI];
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}


void omegaWallFunctionFvPatchScalarField::manipulateMatrix
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


void omegaWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& weights
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    DynamicList<label> constraintCells(weights.size());
    DynamicList<scalar> constraintomega(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& omega
        = dimensionedInternalField();

    label nConstrainedCells = 0;


    forAll(weights, faceI)
    {
        // only set the values if the weights are > tolerance
        if (weights[faceI] > tolerance_)
        {
            nConstrainedCells++;

            label cellI = faceCells[faceI];

            constraintCells.append(cellI);
            constraintomega.append(omega[cellI]);
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
        scalarField(constraintomega.xfer())
    );

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void omegaWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchField<scalar>::write(os);
    writeLocalEntries(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    omegaWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
