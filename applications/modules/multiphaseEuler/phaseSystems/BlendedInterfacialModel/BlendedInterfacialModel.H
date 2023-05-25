/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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

Class
    Foam::BlendedInterfacialModel

Description

SourceFiles
    BlendedInterfacialModel.C

\*---------------------------------------------------------------------------*/

#ifndef BlendedInterfacialModel_H
#define BlendedInterfacialModel_H

#include "blendingMethod.H"
#include "phaseInterface.H"
#include "HashPtrTable.H"
#include "hashedWordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class BlendedInterfacialModel Declaration
\*---------------------------------------------------------------------------*/

template<class ModelType>
class BlendedInterfacialModel
:
    public regIOobject
{
    // Private Data

        //- The interface
        const phaseInterface interface_;

        //- Blending method
        autoPtr<blendingMethod> blending_;

        //- Model for general configurations
        autoPtr<ModelType> modelGeneral_;

        //- Model for phase 1 dispersed in phase 2
        autoPtr<ModelType> model1DispersedIn2_;

        //- Model for phase 2 dispersed in phase 1
        autoPtr<ModelType> model2DispersedIn1_;

        //- Model for phase 1 segregated with phase 2
        autoPtr<ModelType> model1SegregatedWith2_;

        //- Models for general configurations displaced by a third phase
        PtrList<ModelType> modelsGeneralDisplaced_;

        //- Models for phase 1 dispersed in phase 2 displaced by a third phase
        PtrList<ModelType> models1DispersedIn2Displaced_;

        //- Models for phase 2 dispersed in phase 1 displaced by a third phase
        PtrList<ModelType> models2DispersedIn1Displaced_;

        //- Models for phase 1 segregated with phase 2 displaced by a third
        //  phase
        PtrList<ModelType> models1SegregatedWith2Displaced_;

        //- Time index of last check
        mutable label checkTimeIndex_;


    // Private Member Functions

        //- Check compatibility of the available models and the blending method
        void check() const;

        //- Calculate the blending coefficients
        template<template<class> class PatchField, class GeoMesh>
        void calculateBlendingCoeffs
        (
            const UPtrList<const volScalarField>& alphas,
            tmp<GeometricField<scalar, PatchField, GeoMesh>>& fG,
            tmp<GeometricField<scalar, PatchField, GeoMesh>>& f1D2,
            tmp<GeometricField<scalar, PatchField, GeoMesh>>& f2D1,
            tmp<GeometricField<scalar, PatchField, GeoMesh>>& fS,
            PtrList<GeometricField<scalar, PatchField, GeoMesh>>& fGD,
            PtrList<GeometricField<scalar, PatchField, GeoMesh>>& f1D2D,
            PtrList<GeometricField<scalar, PatchField, GeoMesh>>& f2D1D,
            PtrList<GeometricField<scalar, PatchField, GeoMesh>>& fSD,
            const bool subtract
        ) const;

        //- Correct coeff/value on fixed flux boundary conditions
        template<class Type, template<class> class PatchField, class GeoMesh>
        void correctFixedFluxBCs
        (
            GeometricField<Type, PatchField, GeoMesh>& field
        ) const;

        //- Write out a set or surface that visualises the various models'
        //  utilisation of the blending space
        void postProcessBlendingCoefficients(const word& format) const;


protected:

    // Protected Member Functions

        //- Return a blended field
        template
        <
            class Type,
            template<class> class PatchField,
            class GeoMesh,
            class ... Args
        >
        tmp<GeometricField<Type, PatchField, GeoMesh>> evaluate
        (
            tmp<GeometricField<Type, PatchField, GeoMesh>>
            (ModelType::*method)(Args ...) const,
            const word& name,
            const dimensionSet& dims,
            const bool subtract,
            Args ... args
        ) const;

        //- Return a table of blended fields
        template
        <
            class Type,
            template<class> class PatchField,
            class GeoMesh,
            class ... Args
        >
        HashPtrTable<GeometricField<Type, PatchField, GeoMesh>> evaluate
        (
            HashPtrTable<GeometricField<Type, PatchField, GeoMesh>>
            (ModelType::*method)(Args ...) const,
            const word& name,
            const dimensionSet& dims,
            const bool subtract,
            Args ... args
        ) const;

        //- Return a bool combined (or) from all models
        template<class ... Args>
        bool evaluate
        (
            bool (ModelType::*method)(Args ...) const,
            Args ... args
        ) const;

        //- Return a hashed word list combined from all models
        template<class ... Args>
        hashedWordList evaluate
        (
            const hashedWordList& (ModelType::*method)(Args ...) const,
            Args ... args
        ) const;


public:

    //- Runtime type information
    TypeName("BlendedInterfacialModel");


    // Constructors

        //- Construct from a dictionary and an interface
        BlendedInterfacialModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );

        //- Disallow default bitwise copy construction
        BlendedInterfacialModel
        (
            const BlendedInterfacialModel<ModelType>&
        ) = delete;


    //- Destructor
    ~BlendedInterfacialModel();


    // Member Functions

        //- Access the interface
        const phaseInterface& interface() const;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const BlendedInterfacialModel<ModelType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineBlendedInterfacialModelTypeNameAndDebug(ModelType, DebugSwitch)  \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        BlendedInterfacialModel<ModelType>,                                    \
        (                                                                      \
            word(BlendedInterfacialModel<ModelType>::typeName_()) + "<"        \
          + ModelType::typeName + ">"                                          \
        ).c_str(),                                                             \
        DebugSwitch                                                            \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BlendedInterfacialModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
