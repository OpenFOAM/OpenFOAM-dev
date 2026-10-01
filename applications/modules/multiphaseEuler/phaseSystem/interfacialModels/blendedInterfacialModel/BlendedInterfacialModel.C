/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2026 OpenFOAM Foundation
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
#include "generateInterfacialModels.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ModelType>
template<class Type, class GeoMesh, class Method, class ... Args>
Foam::tmp<Foam::GeometricField<Type, GeoMesh>>
Foam::BlendedInterfacialModel<ModelType>::evaluateValueOrField
(
    const Method method,
    const word& name,
    const dimensionSet& dims,
    const Args& ... args
) const
{
    check(ModelType::typeName, models_);

    typedef GeometricField<scalar, GeoMesh> scalarGeoField;
    typedef GeometricField<Type, GeoMesh> typeGeoField;

    const phaseSystem& fluid = interface_.fluid();
    const label nPhases = fluid.phases().size();

    // Get the blending coefficients
    TmpSet<scalarGeoField> fs(nPhases);
    tmp<scalarGeoField>& fG = fs.general;
    tmp<scalarGeoField>& f1D2 = fs.oneDispersedInTwo;
    tmp<scalarGeoField>& f2D1 = fs.twoDispersedInOne;
    tmp<scalarGeoField>& fS = fs.segregated;
    PtrList<scalarGeoField>& fGD = fs.generalDisplaced;
    PtrList<scalarGeoField>& f1D2D = fs.oneDispersedInTwoDisplaced;
    PtrList<scalarGeoField>& f2D1D = fs.twoDispersedInOneDisplaced;
    PtrList<scalarGeoField>& fSD = fs.segregatedDisplaced;
    this->calculateBlendingCoeffs
    (
        fluid.phases().PtrList<phaseModel>::convert<const volScalarField>(),
        models_,
        fs
    );

    // Construct the result
    tmp<typeGeoField> x = typeGeoField::New
    (
        word
        (
            ModelType::typeName, ':',
            IOobject::groupName(name, interface_.name())
        ),
        interface_.mesh(),
        dimensioned<Type>(dims, Zero)
    );

    // Add the model contributions to the result
    if (fG.valid() && models_.general.valid())
    {
        x.ref() += fG*(models_.general().*method)(args ...);
    }
    if (f1D2.valid() && models_.oneDispersedInTwo.valid())
    {
        x.ref() += f1D2*(models_.oneDispersedInTwo().*method)(args ...);
    }
    if (f2D1.valid() && models_.twoDispersedInOne.valid())
    {
        x.ref() += f2D1*(models_.twoDispersedInOne().*method)(args ...);
    }
    if (fS.valid() && models_.segregated.valid())
    {
        x.ref() += fS*(models_.segregated().*method)(args ...);
    }

    forAll(interface_.fluid().phases(), phasei)
    {
        if (fGD.set(phasei) && models_.generalDisplaced.set(phasei))
        {
            x.ref() +=
                fGD[phasei]
               *(models_.generalDisplaced[phasei].*method)(args ...);
        }
        if (f1D2D.set(phasei) && models_.oneDispersedInTwoDisplaced.set(phasei))
        {
            x.ref() +=
                f1D2D[phasei]
               *(models_.oneDispersedInTwoDisplaced[phasei].*method)(args ...);
        }
        if (f2D1D.set(phasei) && models_.twoDispersedInOneDisplaced.set(phasei))
        {
            x.ref() +=
                f2D1D[phasei]
               *(models_.twoDispersedInOneDisplaced[phasei].*method)(args ...);
        }
        if (fSD.set(phasei) && models_.segregatedDisplaced.set(phasei))
        {
            x.ref() +=
                fSD[phasei]
               *(models_.segregatedDisplaced[phasei].*method)(args ...);
        }
    }

    return x;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class ModelType>
template<class Type, class GeoMesh, class ... Args>
Foam::tmp<Foam::GeometricField<Type, GeoMesh>>
Foam::BlendedInterfacialModel<ModelType>::evaluate
(
    const dimensioned<Type>& (ModelType::*method)(Args ...) const,
    const word& name,
    const dimensionSet& dims,
    const Args& ... args
) const
{
    return evaluateValueOrField<Type, GeoMesh>(method, name, dims, args ...);
}


template<class ModelType>
template<class Type, class GeoMesh, class ... Args>
Foam::tmp<Foam::GeometricField<Type, GeoMesh>>
Foam::BlendedInterfacialModel<ModelType>::evaluate
(
    tmp<GeometricField<Type, GeoMesh>> (ModelType::*method)(Args ...) const,
    const word& name,
    const dimensionSet& dims,
    const Args& ... args
) const
{
    return evaluateValueOrField<Type, GeoMesh>(method, name, dims, args ...);
}


template<class ModelType>
template<class Type, class GeoMesh, class ... Args>
Foam::HashPtrTable<Foam::GeometricField<Type, GeoMesh>>
Foam::BlendedInterfacialModel<ModelType>::evaluate
(
    HashPtrTable<GeometricField<Type, GeoMesh>>
    (ModelType::*method)(Args ...) const,
    const word& name,
    const dimensionSet& dims,
    const Args& ... args
) const
{
    check(ModelType::typeName, models_);

    typedef GeometricField<scalar, GeoMesh> scalarGeoField;
    typedef GeometricField<Type, GeoMesh> typeGeoField;

    const phaseSystem& fluid = interface_.fluid();
    const label nPhases = fluid.phases().size();

    // Get the blending coefficients
    TmpSet<scalarGeoField> fs(nPhases);
    tmp<scalarGeoField>& fG = fs.general;
    tmp<scalarGeoField>& f1D2 = fs.oneDispersedInTwo;
    tmp<scalarGeoField>& f2D1 = fs.twoDispersedInOne;
    tmp<scalarGeoField>& fS = fs.segregated;
    PtrList<scalarGeoField>& fGD = fs.generalDisplaced;
    PtrList<scalarGeoField>& f1D2D = fs.oneDispersedInTwoDisplaced;
    PtrList<scalarGeoField>& f2D1D = fs.twoDispersedInOneDisplaced;
    PtrList<scalarGeoField>& fSD = fs.segregatedDisplaced;
    calculateBlendingCoeffs
    (
        fluid.phases().PtrList<phaseModel>::convert<const volScalarField>(),
        models_,
        fs
    );

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
                        word
                        (
                            ModelType::typeName, ':',
                            IOobject::groupName
                            (
                                IOobject::groupName(name, dxIter.key()),
                                interface_.name()
                            )
                        ),
                        f**dxIter()
                    ).ptr()
                );
            }
        }
    };

    if (fG.valid() && models_.general.valid())
    {
        addToXs(fG, (models_.general().*method)(args ...));
    }
    if (f1D2.valid() && models_.oneDispersedInTwo.valid())
    {
        addToXs(f1D2, (models_.oneDispersedInTwo().*method)(args ...));
    }
    if (f2D1.valid() && models_.twoDispersedInOne.valid())
    {
        addToXs(f2D1, (models_.twoDispersedInOne().*method)(args ...));
    }
    if (fS.valid() && models_.segregated.valid())
    {
        addToXs(fS, (models_.segregated().*method)(args ...));
    }

    forAll(interface_.fluid().phases(), phasei)
    {
        if (fGD.set(phasei) && models_.generalDisplaced.set(phasei))
        {
            addToXs
            (
                fGD[phasei],
                (models_.generalDisplaced[phasei].*method)(args ...)
            );
        }
        if (f1D2D.set(phasei) && models_.oneDispersedInTwoDisplaced.set(phasei))
        {
            addToXs
            (
                f1D2D[phasei],
                (models_.oneDispersedInTwoDisplaced[phasei].*method)(args ...)
            );
        }
        if (f2D1D.set(phasei) && models_.twoDispersedInOneDisplaced.set(phasei))
        {
            addToXs
            (
                f2D1D[phasei],
                (models_.twoDispersedInOneDisplaced[phasei].*method)(args ...)
            );
        }
        if (fSD.set(phasei) && models_.segregatedDisplaced.set(phasei))
        {
            addToXs
            (
                fSD[phasei],
                (models_.segregatedDisplaced[phasei].*method)(args ...)
            );
        }
    }

    return xs;
}


template<class ModelType>
template<class ... Args>
bool Foam::BlendedInterfacialModel<ModelType>::evaluate
(
    bool (ModelType::*method)(Args ...) const,
    const Args& ... args
) const
{
    check(ModelType::typeName, models_);

    bool result = false;

    if (models_.general.valid())
    {
        result = result || (models_.general().*method)(args ...);
    }
    if (models_.oneDispersedInTwo.valid())
    {
        result = result || (models_.oneDispersedInTwo().*method)(args ...);
    }
    if (models_.twoDispersedInOne.valid())
    {
        result = result || (models_.twoDispersedInOne().*method)(args ...);
    }
    if (models_.segregated.valid())
    {
        result = result || (models_.segregated().*method)(args ...);
    }

    forAll(interface_.fluid().phases(), phasei)
    {
        if (models_.generalDisplaced.set(phasei))
        {
            result =
                result
             || (models_.generalDisplaced[phasei].*method)(args ...);
        }
        if (models_.oneDispersedInTwoDisplaced.set(phasei))
        {
            result =
                result
             || (models_.oneDispersedInTwoDisplaced[phasei].*method)(args ...);
        }
        if (models_.twoDispersedInOneDisplaced.set(phasei))
        {
            result =
                result
             || (models_.twoDispersedInOneDisplaced[phasei].*method)(args ...);
        }
        if (models_.segregatedDisplaced.set(phasei))
        {
            result =
                result
              || (models_.segregatedDisplaced[phasei].*method)(args ...);
        }
    }

    return result;
}


template<class ModelType>
template<class ... Args>
Foam::hashedWordList Foam::BlendedInterfacialModel<ModelType>::evaluate
(
    const hashedWordList& (ModelType::*method)(Args ...) const,
    const Args& ... args
) const
{
    check(ModelType::typeName, models_);

    wordList result;

    if (models_.general.valid())
    {
        result.append((models_.general().*method)(args ...));
    }
    if (models_.oneDispersedInTwo.valid())
    {
        result.append((models_.oneDispersedInTwo().*method)(args ...));
    }
    if (models_.twoDispersedInOne.valid())
    {
        result.append((models_.twoDispersedInOne().*method)(args ...));
    }
    if (models_.segregated.valid())
    {
        result.append((models_.segregated().*method)(args ...));
    }

    forAll(interface_.fluid().phases(), phasei)
    {
        if (models_.generalDisplaced.set(phasei))
        {
            result.append
            (
                (models_.generalDisplaced[phasei].*method)(args ...)
            );
        }
        if (models_.oneDispersedInTwoDisplaced.set(phasei))
        {
            result.append
            (
                (models_.oneDispersedInTwoDisplaced[phasei].*method)(args ...)
            );
        }
        if (models_.twoDispersedInOneDisplaced.set(phasei))
        {
            result.append
            (
                (models_.twoDispersedInOneDisplaced[phasei].*method)(args ...)
            );
        }
        if (models_.segregatedDisplaced.set(phasei))
        {
            result.append
            (
                 (models_.segregatedDisplaced[phasei].*method)(args ...)
            );
        }
    }

    return hashedWordList(move(result));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ModelType>
template<class ... Args>
Foam::BlendedInterfacialModel<ModelType>::BlendedInterfacialModel
(
    const UPtrList<const entry>& entries,
    const phaseInterface& interface,
    const entry& blendingEntry,
    const Args& ... args
)
:
    blendedInterfacialModel
    (
        typeName,
        ModelType::typeName,
        interface,
        blendingEntry
    ),
    models_(interface.fluid().phases().size())
{
    // Construct the models
    PtrList<phaseInterface> interfaces;
    PtrList<ModelType> models;
    generateInterfacialModels
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
        interfaces,
        models,
        interface.fluid(),
        entries,
        wordHashSet({"blending"}),
        interface,
        args ...
    );

    // Unpack the interface and model lists to populate the models used for the
    // different parts of the blending space
    forAll(interfaces, i)
    {
        const phaseInterface& interface = interfaces[i];

        autoPtr<ModelType>* modelPtrPtr;
        PtrList<ModelType>* modelPtrsPtr;

        if (isA<dispersedPhaseInterface>(interface))
        {
            const phaseModel& dispersed =
                refCast<const dispersedPhaseInterface>(interface).dispersed();

            modelPtrPtr =
                interface_.index(dispersed) == 0
              ? &models_.oneDispersedInTwo
              : &models_.twoDispersedInOne;
            modelPtrsPtr =
                interface_.index(dispersed) == 0
              ? &models_.oneDispersedInTwoDisplaced
              : &models_.twoDispersedInOneDisplaced;
        }
        else if (isA<segregatedPhaseInterface>(interface))
        {
            modelPtrPtr = &models_.segregated;
            modelPtrsPtr = &models_.segregatedDisplaced;
        }
        else
        {
            modelPtrPtr = &models_.general;
            modelPtrsPtr = &models_.generalDisplaced;
        }

        if (!isA<displacedPhaseInterface>(interface))
        {
            *modelPtrPtr = models.set(i, nullptr);
        }
        else
        {
            const phaseModel& displacing =
                refCast<const displacedPhaseInterface>(interface).displacing();

            modelPtrsPtr->set(displacing.index(), models.set(i, nullptr));
        }
    }

    // Write the coefficients to disk
    if (debug)
    {
        this->writeBlendingCoefficients
        (
            ModelType::typeName,
            models_
        );
    }

    // Write a graph or surface of the blending space
    this->postProcessBlendingCoefficients
    (
        ModelType::typeName,
        models_,
        blendingEntry
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ModelType>
Foam::BlendedInterfacialModel<ModelType>::~BlendedInterfacialModel()
{}


// ************************************************************************* //
