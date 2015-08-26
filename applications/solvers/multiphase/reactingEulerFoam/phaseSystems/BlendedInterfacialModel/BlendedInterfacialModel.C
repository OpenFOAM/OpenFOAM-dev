/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 OpenFOAM Foundation
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

#include "BlendedInterfacialModel.H"
#include "fixedValueFvsPatchFields.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ModelType>
template<class GeometricField>
void Foam::BlendedInterfacialModel<ModelType>::correctFixedFluxBCs
(
    GeometricField& field
) const
{
    forAll(phase1_.phi()->boundaryField(), patchI)
    {
        if
        (
            isA<fixedValueFvsPatchScalarField>
            (
                phase1_.phi()->boundaryField()[patchI]
            )
        )
        {
            field.boundaryField()[patchI]
              = pTraits<typename GeometricField::value_type>::zero;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ModelType>
Foam::BlendedInterfacialModel<ModelType>::BlendedInterfacialModel
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const blendingMethod& blending,
    autoPtr<ModelType> model,
    autoPtr<ModelType> model1In2,
    autoPtr<ModelType> model2In1,
    const bool correctFixedFluxBCs
)
:
    phase1_(phase1),
    phase2_(phase2),
    blending_(blending),
    model_(model),
    model1In2_(model1In2),
    model2In1_(model2In1),
    correctFixedFluxBCs_(correctFixedFluxBCs)
{}


template<class ModelType>
Foam::BlendedInterfacialModel<ModelType>::BlendedInterfacialModel
(
    const phasePair::dictTable& modelTable,
    const blendingMethod& blending,
    const phasePair& pair,
    const orderedPhasePair& pair1In2,
    const orderedPhasePair& pair2In1,
    const bool correctFixedFluxBCs
)
:
    phase1_(pair.phase1()),
    phase2_(pair.phase2()),
    blending_(blending),
    correctFixedFluxBCs_(correctFixedFluxBCs)
{
    if (modelTable.found(pair))
    {
        model_.set
        (
            ModelType::New
            (
                modelTable[pair],
                pair
            ).ptr()
        );
    }

    if (modelTable.found(pair1In2))
    {
        model1In2_.set
        (
            ModelType::New
            (
                modelTable[pair1In2],
                pair1In2
            ).ptr()
        );
    }

    if (modelTable.found(pair2In1))
    {
        model2In1_.set
        (
            ModelType::New
            (
                modelTable[pair2In1],
                pair2In1
            ).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ModelType>
Foam::BlendedInterfacialModel<ModelType>::~BlendedInterfacialModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ModelType>
bool Foam::BlendedInterfacialModel<ModelType>::hasModel
(
    const class phaseModel& phase
) const
{
    return
       &phase == &(phase1_)
      ? model1In2_.valid()
      : model2In1_.valid();
}


template<class ModelType>
const ModelType& Foam::BlendedInterfacialModel<ModelType>::model
(
    const class phaseModel& phase
) const
{
    return &phase == &(phase1_) ? model1In2_ : model2In1_;
}


template<class ModelType>
Foam::tmp<Foam::volScalarField>
Foam::BlendedInterfacialModel<ModelType>::K() const
{
    tmp<volScalarField> f1, f2;

    if (model_.valid() || model1In2_.valid())
    {
        f1 = blending_.f1(phase1_, phase2_);
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 = blending_.f2(phase1_, phase2_);
    }

    tmp<volScalarField> x
    (
        new volScalarField
        (
            IOobject
            (
                ModelType::typeName + ":K",
                phase1_.mesh().time().timeName(),
                phase1_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase1_.mesh(),
            dimensionedScalar("zero", ModelType::dimK, 0)
        )
    );

    if (model_.valid())
    {
        x() += model_->K()*(scalar(1) - f1() - f2());
    }
    if (model1In2_.valid())
    {
        x() += model1In2_->K()*f1;
    }
    if (model2In1_.valid())
    {
        x() += model2In1_->K()*f2;
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        correctFixedFluxBCs(x());
    }

    return x;
}


template<class ModelType>
Foam::tmp<Foam::volScalarField>
Foam::BlendedInterfacialModel<ModelType>::K(const scalar residualAlpha) const
{
    tmp<volScalarField> f1, f2;

    if (model_.valid() || model1In2_.valid())
    {
        f1 = blending_.f1(phase1_, phase2_);
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 = blending_.f2(phase1_, phase2_);
    }

    tmp<volScalarField> x
    (
        new volScalarField
        (
            IOobject
            (
                ModelType::typeName + ":K",
                phase1_.mesh().time().timeName(),
                phase1_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase1_.mesh(),
            dimensionedScalar("zero", ModelType::dimK, 0)
        )
    );

    if (model_.valid())
    {
        x() += model_->K(residualAlpha)*(scalar(1) - f1() - f2());
    }
    if (model1In2_.valid())
    {
        x() += model1In2_->K(residualAlpha)*f1;
    }
    if (model2In1_.valid())
    {
        x() += model2In1_->K(residualAlpha)*f2;
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        correctFixedFluxBCs(x());
    }

    return x;
}


template<class ModelType>
Foam::tmp<Foam::surfaceScalarField>
Foam::BlendedInterfacialModel<ModelType>::Kf() const
{
    tmp<surfaceScalarField> f1, f2;

    if (model_.valid() || model1In2_.valid())
    {
        f1 = fvc::interpolate
        (
            blending_.f1(phase1_, phase2_)
        );
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 = fvc::interpolate
        (
            blending_.f2(phase1_, phase2_)
        );
    }

    tmp<surfaceScalarField> x
    (
        new surfaceScalarField
        (
            IOobject
            (
                ModelType::typeName + ":Kf",
                phase1_.mesh().time().timeName(),
                phase1_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase1_.mesh(),
            dimensionedScalar("zero", ModelType::dimK, 0)
        )
    );

    if (model_.valid())
    {
        x() += model_->Kf()*(scalar(1) - f1() - f2());
    }

    if (model1In2_.valid())
    {
        x() += model1In2_->Kf()*f1;
    }

    if (model2In1_.valid())
    {
        x() += model2In1_->Kf()*f2;
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        correctFixedFluxBCs(x());
    }

    return x;
}


template<class ModelType>
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::BlendedInterfacialModel<ModelType>::F() const
{
    tmp<volScalarField> f1, f2;

    if (model_.valid() || model1In2_.valid())
    {
        f1 = blending_.f1(phase1_, phase2_);
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 = blending_.f2(phase1_, phase2_);
    }

    tmp<GeometricField<Type, fvPatchField, volMesh> > x
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                ModelType::typeName + ":F",
                phase1_.mesh().time().timeName(),
                phase1_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase1_.mesh(),
            dimensioned<Type>("zero", ModelType::dimF, pTraits<Type>::zero)
        )
    );

    if (model_.valid())
    {
        x() += model_->F()*(scalar(1) - f1() - f2());
    }

    if (model1In2_.valid())
    {
        x() += model1In2_->F()*f1;
    }

    if (model2In1_.valid())
    {
        x() -= model2In1_->F()*f2; // note : subtraction
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        correctFixedFluxBCs(x());
    }

    return x;
}


template<class ModelType>
Foam::tmp<Foam::surfaceScalarField>
Foam::BlendedInterfacialModel<ModelType>::Ff() const
{
    tmp<surfaceScalarField> f1, f2;

    if (model_.valid() || model1In2_.valid())
    {
        f1 = fvc::interpolate
        (
            blending_.f1(phase1_, phase2_)
        );
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 = fvc::interpolate
        (
            blending_.f2(phase1_, phase2_)
        );
    }

    tmp<surfaceScalarField> x
    (
        new surfaceScalarField
        (
            IOobject
            (
                ModelType::typeName + ":Ff",
                phase1_.mesh().time().timeName(),
                phase1_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase1_.mesh(),
            dimensionedScalar("zero", ModelType::dimF*dimArea, 0)
        )
    );

    if (model_.valid())
    {
        x() += model_->Ff()*(scalar(1) - f1() - f2());
    }

    if (model1In2_.valid())
    {
        x() += model1In2_->Ff()*f1;
    }

    if (model2In1_.valid())
    {
        x() -= model2In1_->Ff()*f2; // note : subtraction
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        correctFixedFluxBCs(x());
    }

    return x;
}


template<class ModelType>
Foam::tmp<Foam::volScalarField>
Foam::BlendedInterfacialModel<ModelType>::D() const
{
    tmp<volScalarField> f1, f2;

    if (model_.valid() || model1In2_.valid())
    {
        f1 = blending_.f1(phase1_, phase2_);
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 = blending_.f2(phase1_, phase2_);
    }

    tmp<volScalarField> x
    (
        new volScalarField
        (
            IOobject
            (
                ModelType::typeName + ":D",
                phase1_.mesh().time().timeName(),
                phase1_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase1_.mesh(),
            dimensionedScalar("zero", ModelType::dimD, 0)
        )
    );

    if (model_.valid())
    {
        x() += model_->D()*(scalar(1) - f1() - f2());
    }
    if (model1In2_.valid())
    {
        x() += model1In2_->D()*f1;
    }
    if (model2In1_.valid())
    {
        x() += model2In1_->D()*f2;
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        correctFixedFluxBCs(x());
    }

    return x;
}


// ************************************************************************* //
