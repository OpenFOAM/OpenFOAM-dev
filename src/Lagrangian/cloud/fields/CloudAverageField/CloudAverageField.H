/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::CloudAverageField

Description
    A cloud field that is averaged and then interpolated back to the cloud. Uses
    CloudDerivedField to provide flexible access and caching.

See also
    Foam::CloudDerivedField

SourceFiles
    CloudAverageField.C

\*---------------------------------------------------------------------------*/

#ifndef CloudAverageField_H
#define CloudAverageField_H

#include "CloudStateField.H"
#include "CloudDerivedField.H"
#include "LagrangianAverage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class volMesh;

/*---------------------------------------------------------------------------*\
                        Class CloudAverageField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CloudAverageField
:
    public CloudDerivedField<Type>
{
    // Private Static Data

        //- Pointer to remove method for constant weight sum
        static constexpr void (LagrangianAverage<Type>::*removeWeightSum_)
        (
            const LagrangianSubSubField<Type>&
        ) = &LagrangianAverage<Type>::remove;

        //- Pointer to remove method for computed weight sum
        static constexpr void (LagrangianAverage<Type>::*removeNoWeightSum_)
        (
            const LagrangianSubSubField<scalar>&,
            const LagrangianSubSubField<Type>&
        ) = &LagrangianAverage<Type>::remove;

        //- Pointer to add method for constant weight sum
        static constexpr void (LagrangianAverage<Type>::*addWeightSum_)
        (
            const LagrangianSubSubField<Type>&,
            const bool
        ) = &LagrangianAverage<Type>::add;

        //- Pointer to add method for computed weight sum
        static constexpr void (LagrangianAverage<Type>::*addNoWeightSum_)
        (
            const LagrangianSubSubField<scalar>&,
            const LagrangianSubSubField<Type>&,
            const bool
        ) = &LagrangianAverage<Type>::add;

        //- Pointer to correct method for constant weight sum
        static constexpr void (LagrangianAverage<Type>::*correctWeightSum_)
        (
            const LagrangianSubSubField<Type>&,
            const bool
        ) = &LagrangianAverage<Type>::correct;

        //- Pointer to correct method for computed weight sum
        static constexpr void (LagrangianAverage<Type>::*correctNoWeightSum_)
        (
            const LagrangianSubSubField<scalar>&,
            const LagrangianSubSubField<Type>&,
            const bool
        ) = &LagrangianAverage<Type>::correct;


    // Private Data

        //- The name of the average
        const word name_;

        //- The cell weight sum
        const tmp<DimensionedField<scalar, volMesh>> tcellWeightSum_;

        //- Reference to the weight state field
        const CloudStateField<scalar>& weightState_;

        //- Reference to the weight derived field
        const CloudDerivedField<scalar>& weightDerived_;

        //- Reference to the state field
        const CloudStateField<Type>& weightPsiOrPsiState_;

        //- Reference to the derived field
        const CloudDerivedField<Type>& weightPsiOrPsiDerived_;

        //- Averaging engine
        mutable autoPtr<LagrangianAverage<Type>> psiAveragePtr_;

        //- Current state of the average
        mutable enum class psiAverageState
        {
            removed,
            complete,
            cached
        } psiAverageState_;


    // Private Member Functions

        //- Run an operation on the averaging engine
        template<class MethodWeightSum, class MethodNoWeightSum, class ... Args>
        void op
        (
            MethodWeightSum,
            MethodNoWeightSum,
            const LagrangianSubMesh&,
            const Args& ... args
        ) const;

        //- Interpolate and return the result
        tmp<LagrangianSubField<Type>> interpolate
        (
            const LagrangianModelRef&,
            const LagrangianSubMesh&
        ) const;


public:

    // Constructors

        //- Construct from a name, a cell weight sum and a derived field
        CloudAverageField
        (
            const word& name,
            const DimensionedField<scalar, volMesh>& cellWeightSum,
            const CloudDerivedField<Type>& weightPsi
        );


    // Member Functions

        //- Remove this sub-mesh from the average
        void remove(const LagrangianSubMesh& subMesh);

        //- Add this sub-mesh to the average
        void add(const LagrangianSubMesh& subMesh, const bool cache);

        //- Correct the average with values from the sub-mesh
        void correct(const LagrangianSubMesh& subMesh, const bool cache);

        //- Reset
        void reset();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CloudAverageField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
