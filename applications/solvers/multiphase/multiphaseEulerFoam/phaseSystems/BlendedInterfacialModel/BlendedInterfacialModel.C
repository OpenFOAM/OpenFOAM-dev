/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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
#include "phaseSystem.H"
#include "dispersedDisplacedPhaseInterface.H"
#include "segregatedDisplacedPhaseInterface.H"
#include "fixedValueFvsPatchFields.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendedInterfacialModel2
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

} // End namespace blendedInterfacialModel2
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ModelType>
template<template<class> class PatchField, class GeoMesh>
void Foam::BlendedInterfacialModel<ModelType>::calculateBlendingCoeffs
(
    tmp<GeometricField<scalar, PatchField, GeoMesh>>& f1D2,
    tmp<GeometricField<scalar, PatchField, GeoMesh>>& f2D1,
    const bool subtract
) const
{
    typedef GeometricField<scalar, PatchField, GeoMesh> scalarGeoField;

    if ((modelGeneral_.valid() || model1SegregatedWith2_.valid()) && subtract)
    {
        FatalErrorInFunction
            << "Cannot treat an interfacial model with no distinction between "
            << "continuous and dispersed phases as signed"
            << exit(FatalError);
    }

    if (model1SegregatedWith2_.valid() || model1DispersedIn2_.valid())
    {
        f1D2 =
            blendedInterfacialModel2::interpolate<scalarGeoField>
            (
                blending_->f1(interface_.phase1(), interface_.phase2())
            );
    }

    if (model1SegregatedWith2_.valid() || model2DispersedIn1_.valid())
    {
        f2D1 =
            (subtract ? -1 : +1)
           *blendedInterfacialModel2::interpolate<scalarGeoField>
            (
                blending_->f2(interface_.phase1(), interface_.phase2())
            );
    }
}


template<class ModelType>
template<template<class> class PatchField, class GeoMesh>
void Foam::BlendedInterfacialModel<ModelType>::calculateBlendingCoeffs
(
    tmp<GeometricField<scalar, PatchField, GeoMesh>>& fG,
    tmp<GeometricField<scalar, PatchField, GeoMesh>>& f1D2,
    tmp<GeometricField<scalar, PatchField, GeoMesh>>& f2D1,
    tmp<GeometricField<scalar, PatchField, GeoMesh>>& fS,
    const bool subtract
) const
{
    typedef GeometricField<scalar, PatchField, GeoMesh> scalarGeoField;

    // Get the dispersed blending coefficients
    calculateBlendingCoeffs(f1D2, f2D1, subtract);

    // Create a segregated blending coefficient if necessary
    if (model1SegregatedWith2_.valid())
    {
        fS =
            scalarGeoField::New
            (
                ModelType::typeName + ":"
              + IOobject::groupName("fS", interface_.name()),
                interface_.mesh(),
                dimensionedScalar(dimless, 1)
            );
        if (model1DispersedIn2_.valid())
        {
            fS.ref() -= f1D2();
        }
        if (model2DispersedIn1_.valid())
        {
            fS.ref() -= f2D1();
        }
    }

    // Create a general blending coefficient if necessary
    if (modelGeneral_.valid())
    {
        fG =
            scalarGeoField::New
            (
                ModelType::typeName + ":"
              + IOobject::groupName("fG", interface_.name()),
                interface_.mesh(),
                dimensionedScalar(dimless, 1)
            );

        if (model1DispersedIn2_.valid())
        {
            fG.ref() -= f1D2();
        }
        if (model2DispersedIn1_.valid())
        {
            fG.ref() -= f2D1();
        }
        if (model1SegregatedWith2_.valid())
        {
            fG.ref() -= fS();
        }
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
                !interface_.phase1().stationary()
             && isA<fixedValueFvsPatchScalarField>
                (
                    interface_.phase1().phi()().boundaryField()[patchi]
                )
            )
         || (
                !interface_.phase2().stationary()
             && isA<fixedValueFvsPatchScalarField>
                (
                    interface_.phase2().phi()().boundaryField()[patchi]
                )
            )
        )
        {
            fieldBf[patchi] = Zero;
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

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

    // Get the blending coefficients
    tmp<scalarGeoField> fG, f1D2, f2D1, fS;
    calculateBlendingCoeffs(fG, f1D2, f2D1, fS, subtract);

    // Construct the result
    tmp<typeGeoField> x =
        typeGeoField::New
        (
            ModelType::typeName + ":"
          + IOobject::groupName(name, interface_.name()),
            interface_.mesh(),
            dimensioned<Type>(dims, Zero)
        );

    // Add the model contributions to the result
    if (modelGeneral_.valid())
    {
        x.ref() += fG*(modelGeneral_().*method)(args ...);
    }
    if (model1DispersedIn2_.valid())
    {
        x.ref() += f1D2*(model1DispersedIn2_().*method)(args ...);
    }
    if (model2DispersedIn1_.valid())
    {
        x.ref() += f2D1*(model2DispersedIn1_().*method)(args ...);
    }
    if (model1SegregatedWith2_.valid())
    {
        x.ref() += fS*(model1SegregatedWith2_().*method)(args ...);
    }

    // Correct boundary conditions if necessary
    if
    (
        ModelType::correctFixedFluxBCs
     && (
            modelGeneral_.valid()
         || model1DispersedIn2_.valid()
         || model2DispersedIn1_.valid()
         || model1SegregatedWith2_.valid()
        )
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

    // Get the blending coefficients
    tmp<scalarGeoField> fG, f1D2, f2D1, fS;
    calculateBlendingCoeffs(fG, f1D2, f2D1, fS, subtract);

    // Construct the result
    HashPtrTable<typeGeoField> xs;

    // Add the model contributions to the result
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
                            interface_.name()
                        ),
                        f**dxIter()
                    ).ptr()
                );
            }
        }
    };
    if (modelGeneral_.valid())
    {
        addToXs(fG, (modelGeneral_().*method)(args ...));
    }
    if (model1DispersedIn2_.valid())
    {
        addToXs(f1D2, (model1DispersedIn2_().*method)(args ...));
    }
    if (model2DispersedIn1_.valid())
    {
        addToXs(f2D1, (model2DispersedIn1_().*method)(args ...));
    }
    if (model1SegregatedWith2_.valid())
    {
        addToXs(fS, (model1SegregatedWith2_().*method)(args ...));
    }

    // Correct boundary conditions if necessary
    if
    (
        ModelType::correctFixedFluxBCs
     && (
            modelGeneral_.valid()
         || model1DispersedIn2_.valid()
         || model2DispersedIn1_.valid()
         || model1SegregatedWith2_.valid()
        )
    )
    {
        forAllIter(typename HashPtrTable<typeGeoField>, xs, xIter)
        {
            correctFixedFluxBCs(*xIter());
        }
    }

    return xs;
}


template<class ModelType>
template<class ... Args>
bool Foam::BlendedInterfacialModel<ModelType>::evaluate
(
    bool (ModelType::*method)(Args ...) const,
    Args ... args
) const
{
    return
        (
            modelGeneral_.valid()
         && (modelGeneral_().*method)(args ...)
        )
     || (
            model1DispersedIn2_.valid()
         && (model1DispersedIn2_().*method)(args ...)
        )
     || (
            model2DispersedIn1_.valid()
         && (model2DispersedIn1_().*method)(args ...)
        )
     || (
            model1SegregatedWith2_.valid()
         && (model1SegregatedWith2_().*method)(args ...)
        );
}


template<class ModelType>
template<class ... Args>
Foam::hashedWordList Foam::BlendedInterfacialModel<ModelType>::evaluate
(
    const hashedWordList& (ModelType::*method)(Args ...) const,
    Args ... args
) const
{
    wordList result;

    if (modelGeneral_.valid())
    {
        result.append((modelGeneral_().*method)(args ...));
    }
    if (model1DispersedIn2_.valid())
    {
        result.append((model1DispersedIn2_().*method)(args ...));
    }
    if (model2DispersedIn1_.valid())
    {
        result.append((model2DispersedIn1_().*method)(args ...));
    }
    if (model1SegregatedWith2_.valid())
    {
        result.append((model1SegregatedWith2_().*method)(args ...));
    }

    return hashedWordList(move(result));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ModelType>
Foam::BlendedInterfacialModel<ModelType>::BlendedInterfacialModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, interface.name()),
            interface.fluid().mesh().time().timeName(),
            interface.fluid().mesh()
        )
    ),
    interface_(interface)
{
    // Construct blending functions
    const dictionary& blendingDict = interface.fluid().subDict("blending");
    blending_ =
        blendingMethod::New
        (
            ModelType::typeName,
            blendingDict.found(ModelType::typeName)
          ? blendingDict.subDict(ModelType::typeName)
          : blendingDict.subDict("default"),
            interface.fluid().phases().toc()
        );

    // Construct the models
    PtrList<phaseInterface> interfaces;
    PtrList<ModelType> models;
    interface.fluid().generateInterfacialModels
    <
        ModelType,
        dispersedDisplacedPhaseInterface,
        segregatedDisplacedPhaseInterface,
        displacedPhaseInterface,
        dispersedPhaseInterface,
        segregatedPhaseInterface,
        phaseInterface
    >
    (
        dict,
        interface,
        interfaces,
        models
    );

    // Define local set function
    auto set = [&]
    (
        const phaseInterface& interface,
        ModelType* model,
        autoPtr<ModelType>& modelGeneral,
        autoPtr<ModelType>& model1DispersedIn2,
        autoPtr<ModelType>& model2DispersedIn1,
        autoPtr<ModelType>& model1SegregatedWith2
    )
    {
        if (isA<dispersedPhaseInterface>(interface))
        {
            const phaseModel& dispersed =
                refCast<const dispersedPhaseInterface>(interface).dispersed();

            interface_.index(dispersed) == 0
          ? model1DispersedIn2.set(model)
          : model2DispersedIn1.set(model);
        }
        else if (isA<segregatedPhaseInterface>(interface))
        {
            model1SegregatedWith2.set(model);
        }
        else
        {
            modelGeneral.set(model);
        }
    };

    // Unpack the interface and model lists to populate the models used for the
    // different parts of the blending space
    forAll(interfaces, i)
    {
        if (isA<displacedPhaseInterface>(interfaces[i]))
        {
            /*
            const phaseModel& displacing =
                refCast<const displacedPhaseInterface>
                (interfaces[i]).displacing();

            set
            (
                interfaces[i],
                models[i],
                modelGeneralDisplaced_[displacing.index()],
                model1DispersedIn2Displaced_[displacing.index()],
                model2DispersedIn1Displaced_[displacing.index()],
                model1SegregatedWith2Displaced_[displacing.index()]
            );
            */

            NotImplemented;
        }
        else
        {
            set
            (
                interfaces[i],
                models.set(i, nullptr).ptr(),
                modelGeneral_,
                model1DispersedIn2_,
                model2DispersedIn1_,
                model1SegregatedWith2_
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ModelType>
Foam::BlendedInterfacialModel<ModelType>::~BlendedInterfacialModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ModelType>
const Foam::phaseInterface&
Foam::BlendedInterfacialModel<ModelType>::interface() const
{
    return interface_;
}


template<class ModelType>
bool Foam::BlendedInterfacialModel<ModelType>::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
