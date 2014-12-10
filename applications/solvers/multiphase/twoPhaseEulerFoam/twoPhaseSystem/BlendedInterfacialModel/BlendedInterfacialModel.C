/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class modelType>
template<class Type>
void Foam::BlendedInterfacialModel<modelType>::correctFixedFluxBCs
(
    GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    forAll(pair_.phase1().phi().boundaryField(), patchI)
    {
        if
        (
            isA<fixedValueFvsPatchScalarField>
            (
                pair_.phase1().phi().boundaryField()[patchI]
            )
        )
        {
            field.boundaryField()[patchI] = pTraits<Type>::zero;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class modelType>
Foam::BlendedInterfacialModel<modelType>::BlendedInterfacialModel
(
    const phasePair::dictTable& modelTable,
    const blendingMethod& blending,
    const phasePair& pair,
    const orderedPhasePair& pair1In2,
    const orderedPhasePair& pair2In1
)
:
    pair_(pair),
    pair1In2_(pair1In2),
    pair2In1_(pair2In1),
    blending_(blending)
{
    if (modelTable.found(pair_))
    {
        model_.set
        (
            modelType::New
            (
                modelTable[pair_],
                pair_
            ).ptr()
        );
    }

    if (modelTable.found(pair1In2_))
    {
        model1In2_.set
        (
            modelType::New
            (
                modelTable[pair1In2_],
                pair1In2_
            ).ptr()
        );
    }

    if (modelTable.found(pair2In1_))
    {
        model2In1_.set
        (
            modelType::New
            (
                modelTable[pair2In1_],
                pair2In1_
            ).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class modelType>
Foam::BlendedInterfacialModel<modelType>::~BlendedInterfacialModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class modelType>
Foam::tmp<Foam::volScalarField>
Foam::BlendedInterfacialModel<modelType>::K() const
{
    tmp<volScalarField> f1, f2;

    if (model_.valid() || model1In2_.valid())
    {
        f1 = blending_.f1(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 = blending_.f2(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    tmp<volScalarField> x
    (
        new volScalarField
        (
            IOobject
            (
                modelType::typeName + "Coeff",
                pair_.phase1().mesh().time().timeName(),
                pair_.phase1().mesh()
            ),
            pair_.phase1().mesh(),
            dimensionedScalar("zero", modelType::dimK, 0)
        )
    );

    if (model_.valid())
    {
        x() += model_->K()*(f1() - f2());
    }

    if (model1In2_.valid())
    {
        x() += model1In2_->K()*(1 - f1);
    }

    if (model2In1_.valid())
    {
        x() += model2In1_->K()*f2;
    }

    if (model_.valid() || model1In2_.valid() || model2In1_.valid())
    {
        correctFixedFluxBCs(x());
    }

    return x;
}


template<class modelType>
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::BlendedInterfacialModel<modelType>::F() const
{
    tmp<volScalarField> f1, f2;

    if (model_.valid() || model1In2_.valid())
    {
        f1 = blending_.f1(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 = blending_.f2(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    tmp<GeometricField<Type, fvPatchField, volMesh> > x
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                modelType::typeName + "Coeff",
                pair_.phase1().mesh().time().timeName(),
                pair_.phase1().mesh()
            ),
            pair_.phase1().mesh(),
            dimensioned<Type>("zero", modelType::dimF, pTraits<Type>::zero)
        )
    );

    if (model_.valid())
    {
        x() += model_->F()*(f1() - f2());
    }

    if (model1In2_.valid())
    {
        x() += model1In2_->F()*(1 - f1);
    }

    if (model2In1_.valid())
    {
        x() -= model2In1_->F()*f2; // note : subtraction
    }

    if (model_.valid() || model1In2_.valid() || model2In1_.valid())
    {
        correctFixedFluxBCs(x());
    }

    return x;
}


template<class modelType>
const modelType& Foam::BlendedInterfacialModel<modelType>::phaseModel
(
    const class phaseModel& phase
) const
{
    return &phase == &(pair_.phase1()) ? model1In2_ : model2In1_;
}


// ************************************************************************* //
