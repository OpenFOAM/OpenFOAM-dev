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
                         Class BlendedInterfacialModel Declaration
\*---------------------------------------------------------------------------*/

template<class modelType>
class BlendedInterfacialModel
{
private:

    // Private data

        //- Unordered phase pair
        const phasePair& pair_;

        //- Ordered phase pair for dispersed phase 1 in continuous phase 2
        const orderedPhasePair& pair1In2_;

        //- Ordered phase pair for dispersed phase 2 in continuous phase 1
        const orderedPhasePair& pair2In1_;

        //- Model for region with no obvious dispersed phase
        autoPtr<modelType> model_;

        //- Model for dispersed phase 1 in continuous phase 2
        autoPtr<modelType> model1In2_;

        //- Model for dispersed phase 2 in continuous phase 1
        autoPtr<modelType> model2In1_;

        //- Blending model
        const blendingMethod& blending_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlendedInterfacialModel(const BlendedInterfacialModel<modelType>&);

        //- Disallow default bitwise assignment
        void operator=(const BlendedInterfacialModel<modelType>&);

        //- Correct coeff/value on fixed flux boundary conditions
        template<class Type>
        void correctFixedFluxBCs
        (
            GeometricField<Type, fvPatchField, volMesh>& field
        ) const;


public:

    // Constructors

        //- Construct from the model table, dictionary and pairs
        BlendedInterfacialModel
        (
            const phasePair::dictTable& modelTable,
            const blendingMethod& blending,
            const phasePair& pair,
            const orderedPhasePair& pair1In2,
            const orderedPhasePair& pair2In1
        );


    //- Destructor
    ~BlendedInterfacialModel();


    // Member Functions

        //- Return the implicit coefficient
        tmp<volScalarField> K() const;

        //- Return the explicit value
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> > F() const;
        
        //- Return the model for the supplied phase
        const modelType& phaseModel(const phaseModel& phase) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlendedInterfacialModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
