/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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
#include "phasePair.H"
#include "orderedPhasePair.H"
#include "geometricZeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class blendedInterfacialModel Declaration
\*---------------------------------------------------------------------------*/

class blendedInterfacialModel
{
    public:

        //- Convenience function to interpolate blending values. Needs to be
        //  specialised, so can't sit in the templated class.
        template<class GeoField>
        static inline tmp<GeoField> interpolate(tmp<volScalarField> f);
};


/*---------------------------------------------------------------------------*\
                    Class BlendedInterfacialModel Declaration
\*---------------------------------------------------------------------------*/

template<class ModelType>
class BlendedInterfacialModel
:
    public regIOobject
{
    // Private data

        //- Reference to phase 1
        const phaseModel& phase1_;

        //- Reference to phase 2
        const phaseModel& phase2_;

        //- Blending model
        const blendingMethod& blending_;

        //- Model for region with no obvious dispersed phase
        autoPtr<ModelType> model_;

        //- Model for dispersed phase 1 in continuous phase 2
        autoPtr<ModelType> model1In2_;

        //- Model for dispersed phase 2 in continuous phase 1
        autoPtr<ModelType> model2In1_;

        //- If true set coefficients and forces to 0 at fixed-flux BCs
        bool correctFixedFluxBCs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlendedInterfacialModel(const BlendedInterfacialModel<ModelType>&);

        //- Disallow default bitwise assignment
        void operator=(const BlendedInterfacialModel<ModelType>&);

        //- Correct coeff/value on fixed flux boundary conditions
        template<class GeoField>
        void correctFixedFluxBCs(GeoField& field) const;

        //- Return the blended coeff/value
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


public:

    //- Runtime type information
    TypeName("BlendedInterfacialModel");


    // Constructors

        //- Construct from two phases, blending method and three models
        BlendedInterfacialModel
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const blendingMethod& blending,
            autoPtr<ModelType> model,
            autoPtr<ModelType> model1In2,
            autoPtr<ModelType> model2In1,
            const bool correctFixedFluxBCs = true
        );


        //- Construct from the model table, dictionary and pairs
        BlendedInterfacialModel
        (
            const phasePair::dictTable& modelTable,
            const blendingMethod& blending,
            const phasePair& pair,
            const orderedPhasePair& pair1In2,
            const orderedPhasePair& pair2In1,
            const bool correctFixedFluxBCs = true
        );


    //- Destructor
    ~BlendedInterfacialModel();


    // Member Functions

        //- Return true if a model is specified for the supplied phase
        bool hasModel(const phaseModel& phase) const;

        //- Return the model for the supplied phase
        const ModelType& model(const phaseModel& phase) const;

        //- Return the sign of the explicit value for the supplied phase
        scalar sign(const phaseModel& phase) const;

        //- Return the blended force coefficient
        tmp<volScalarField> K() const;

        //- Return the blended force coefficient with a specified residual alpha
        tmp<volScalarField> K(const scalar residualAlpha) const;

        //- Return the face blended force coefficient
        tmp<surfaceScalarField> Kf() const;

        //- Return the blended force
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>> F() const;

        //- Return the face blended force
        tmp<surfaceScalarField> Ff() const;

        //- Return the blended diffusivity
        tmp<volScalarField> D() const;

        //- Return the blended mass transfer rate
        tmp<volScalarField> dmdt() const;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineBlendedInterfacialModelTypeNameAndDebug(ModelType, DebugSwitch)  \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        BlendedInterfacialModel<ModelType>,                                    \
        (                                                                      \
            word(BlendedInterfacialModel<ModelType>::typeName_()) + "<"        \
          + ModelType::typeName_() + ">"                                       \
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
