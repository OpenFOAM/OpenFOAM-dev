/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2020 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendedInterfacialModel
{

template<class GeoField>
inline tmp<GeoField> interpolate(tmp<volScalarField> f);

template<>
inline tmp<Foam::volScalarField> interpolate(tmp<volScalarField> f)
{
    return f;
}

template<>
inline tmp<Foam::surfaceScalarField> interpolate(tmp<volScalarField> f)
{
    return fvc::interpolate(f);
}

} // End namespace blendedInterfacialModel
} // End namespace Foam

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ModelType>
template<template<class> class PatchField, class GeoMesh>
void Foam::BlendedInterfacialModel<ModelType>::calculateBlendingCoeffs
(
    tmp<GeometricField<scalar, PatchField, GeoMesh>>& f1,
    tmp<GeometricField<scalar, PatchField, GeoMesh>>& f2,
    const bool subtract
) const
{
    typedef GeometricField<scalar, PatchField, GeoMesh> scalarGeoField;

    if (model_.valid() && subtract)
    {
        FatalErrorInFunction
            << "Cannot treat an interfacial model with no distinction between "
            << "continuous and dispersed phases as signed"
            << exit(FatalError);
    }

    if (model_.valid() || model1In2_.valid())
    {
        f1 =
            blendedInterfacialModel::interpolate<scalarGeoField>
            (
                blending_.f1(phase1_, phase2_)
            );
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 =
            (subtract ? -1 : +1)
           *blendedInterfacialModel::interpolate<scalarGeoField>
            (
                blending_.f2(phase1_, phase2_)
            );
    }
}


template<class ModelType>
template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::BlendedInterfacialModel<ModelType>::correctFixedFluxBCs
(
    GeometricField<Type, PatchField, GeoMesh>& field
) const
{
    typedef GeometricField<Type, PatchField, GeoMesh> typeGeoField;

    typename typeGeoField::Boundary& fieldBf = field.boundaryFieldRef();

    forAll(fieldBf, patchi)
    {
        if
        (
            (
                !phase1_.stationary()
             && isA<fixedValueFvsPatchScalarField>
                (
                    phase1_.phi()().boundaryField()[patchi]
                )
            )
         || (
                !phase2_.stationary()
             && isA<fixedValueFvsPatchScalarField>
                (
                    phase2_.phi()().boundaryField()[patchi]
                )
            )
        )
        {
            fieldBf[patchi] = Zero;
        }
    }
}


template<class ModelType>
template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh,
    class ... Args
>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::BlendedInterfacialModel<ModelType>::evaluate
(
    tmp<GeometricField<Type, PatchField, GeoMesh>>
    (ModelType::*method)(Args ...) const,
    const word& name,
    const dimensionSet& dims,
    const bool subtract,
    Args ... args
) const
{
    typedef GeometricField<scalar, PatchField, GeoMesh> scalarGeoField;
    typedef GeometricField<Type, PatchField, GeoMesh> typeGeoField;

    tmp<scalarGeoField> f1, f2;
    calculateBlendingCoeffs(f1, f2, subtract);

    tmp<typeGeoField> x =
        typeGeoField::New
        (
            ModelType::typeName + ":"
          + IOobject::groupName(name, phasePair(phase1_, phase2_).name()),
            phase1_.mesh(),
            dimensioned<Type>(dims, Zero)
        );

    if (model_.valid())
    {
        x.ref() += (scalar(1) - f1() - f2())*(model_().*method)(args ...);
    }

    if (model1In2_.valid())
    {
        x.ref() += f1*(model1In2_().*method)(args ...);
    }

    if (model2In1_.valid())
    {
        x.ref() += f2*(model2In1_().*method)(args ...);
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        correctFixedFluxBCs(x.ref());
    }

    return x;
}


template<class ModelType>
template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh,
    class ... Args
>
Foam::HashPtrTable<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::BlendedInterfacialModel<ModelType>::evaluate
(
    HashPtrTable<GeometricField<Type, PatchField, GeoMesh>>
    (ModelType::*method)(Args ...) const,
    const word& name,
    const dimensionSet& dims,
    const bool subtract,
    Args ... args
) const
{
    typedef GeometricField<scalar, PatchField, GeoMesh> scalarGeoField;
    typedef GeometricField<Type, PatchField, GeoMesh> typeGeoField;

    tmp<scalarGeoField> f1, f2;
    calculateBlendingCoeffs(f1, f2, subtract);

    HashPtrTable<typeGeoField> xs;

    auto addToXs = [&]
    (
        const scalarGeoField& f,
        const HashPtrTable<typeGeoField>& dxs
    )
    {
        forAllConstIter(typename HashPtrTable<typeGeoField>, dxs, dxIter)
        {
            if (xs.found(dxIter.key()))
            {
                *xs[dxIter.key()] += f**dxIter();
            }
            else
            {
                xs.insert
                (
                    dxIter.key(),
                    typeGeoField::New
                    (
                        ModelType::typeName + ':'
                      + IOobject::groupName
                        (
                            IOobject::groupName(name, dxIter.key()),
                            phasePair(phase1_, phase2_).name()
                        ),
                        f**dxIter()
                    ).ptr()
                );
            }
        }
    };

    if (model_.valid())
    {
        addToXs(scalar(1) - f1() - f2(), (model_().*method)(args ...));
    }

    if (model1In2_.valid())
    {
        addToXs(f1, (model1In2_().*method)(args ...));
    }

    if (model2In1_.valid())
    {
        addToXs(f2, (model2In1_().*method)(args ...));
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        forAllIter(typename HashPtrTable<typeGeoField>, xs, xIter)
        {
            correctFixedFluxBCs(*xIter());
        }
    }

    return xs;
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
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, phasePair(phase1, phase2).name()),
            phase1.mesh().time().timeName(),
            phase1.mesh()
        )
    ),
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
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh()
        )
    ),
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
Foam::tmp<Foam::volScalarField>
Foam::BlendedInterfacialModel<ModelType>::K() const
{
    tmp<volScalarField> (ModelType::*k)() const = &ModelType::K;

    return evaluate(k, "K", ModelType::dimK, false);
}


template<class ModelType>
Foam::tmp<Foam::volScalarField>
Foam::BlendedInterfacialModel<ModelType>::K(const scalar residualAlpha) const
{
    tmp<volScalarField> (ModelType::*k)(const scalar) const = &ModelType::K;

    return evaluate(k, "K", ModelType::dimK, false, residualAlpha);
}


template<class ModelType>
Foam::tmp<Foam::surfaceScalarField>
Foam::BlendedInterfacialModel<ModelType>::Kf() const
{
    return evaluate(&ModelType::Kf, "Kf", ModelType::dimK, false);
}


template<class ModelType>
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::BlendedInterfacialModel<ModelType>::F() const
{
    return evaluate(&ModelType::F, "F", ModelType::dimF, true);
}


template<class ModelType>
Foam::tmp<Foam::surfaceScalarField>
Foam::BlendedInterfacialModel<ModelType>::Ff() const
{
    return evaluate(&ModelType::Ff, "Ff", ModelType::dimF*dimArea, true);
}


template<class ModelType>
Foam::tmp<Foam::volScalarField>
Foam::BlendedInterfacialModel<ModelType>::D() const
{
    return evaluate(&ModelType::D, "D", ModelType::dimD, false);
}


template<class ModelType>
bool Foam::BlendedInterfacialModel<ModelType>::mixture() const
{
    return
        (model1In2_.valid() && model1In2_->mixture())
     || (model2In1_.valid() && model2In1_->mixture())
     || (model_.valid() && model_->mixture());
}


template<class ModelType>
Foam::tmp<Foam::volScalarField>
Foam::BlendedInterfacialModel<ModelType>::dmdtf() const
{
    return evaluate(&ModelType::dmdtf, "dmdtf", ModelType::dimDmdt, true);
}


template<class ModelType>
Foam::hashedWordList Foam::BlendedInterfacialModel<ModelType>::species() const
{
    wordList species;

    if (model1In2_.valid())
    {
        species.append(model1In2_->species());
    }
    if (model2In1_.valid())
    {
        species.append(model2In1_->species());
    }
    if (model_.valid())
    {
        species.append(model_->species());
    }

    return hashedWordList(move(species));
}


template<class ModelType>
Foam::HashPtrTable<Foam::volScalarField>
Foam::BlendedInterfacialModel<ModelType>::dmidtf() const
{
    return evaluate(&ModelType::dmidtf, "dmidtf", ModelType::dimDmdt, true);
}


template<class ModelType>
bool Foam::BlendedInterfacialModel<ModelType>::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
