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
    Foam::CarrierEqn

Description
    Class to store the coefficients of a term in a carrier equation

SourceFiles
    CarrierEqn.C

\*---------------------------------------------------------------------------*/

#ifndef CarrierEqn_H
#define CarrierEqn_H

#include "CarrierCoeff.H"
#include "LagrangianEqn.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class CarrierEqn Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CarrierEqn
:
    public tmp<CarrierEqn<Type>>::refCount,
    public CarrierEqnBase
{
    // Private Data

        //- Reference to the volume field
        const VolField<Type>& psi_;


    // Private Member Functions

        //- Return the residual
        template<class PsiFieldType>
        tmp<VolInternalField<Type>> residual
        (
            const dimensionSet& dims,
            const word& psiName,
            const PsiFieldType& psi
        ) const;


public:

    // Public Data

        //- Explicit coefficient
        CarrierCoeff<Type, false> Su;

        //- Implicit coefficient
        CarrierCoeff<scalar, true> Sp;


    // Constructors

        //- Construct for a field
        CarrierEqn(const VolField<Type>& psi);

        //- Construct from a field name and a reference to the mesh. The field
        //  reference is not set and the field is assumed to be (dimensionless)
        //  unity. This is for volume sources which are expressed per unit
        //  volume and therefore have no associated field.
        CarrierEqn(const word& psiName, const fvMesh& mesh);

        //- Copy construct
        CarrierEqn(const CarrierEqn<Type>& eqn);

        //- Construct from tmp
        CarrierEqn(const tmp<CarrierEqn<Type>>& tEqn);


    // Member Functions

        //- Clear the coefficient values
        void clear();

        //- Return the residual
        tmp<VolInternalField<Type>> residual(const dimensionSet& dims) const;


    // Member Operators

        //- Addition assignment
        void operator+=(const LagrangianEqn<Type>& lEqn);

        //- Addition assignment
        void operator+=(const tmp<LagrangianEqn<Type>>& tlEqn);

        //- Subtraction assignment
        void operator-=(const LagrangianEqn<Type>& lEqn);

        //- Subtraction assignment
        void operator-=(const tmp<LagrangianEqn<Type>>& tlEqn);
};


#define CARRIER_EQN_SCALAR_OPERATOR(Op)                                        \
                                                                               \
    template<class Type>                                                       \
    tmp<CarrierEqn<Type>> operator Op                                          \
    (                                                                          \
        const CarrierEqn<Type>& eqn,                                           \
        const dimensionedScalar& s                                             \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    tmp<CarrierEqn<Type>> operator Op                                          \
    (                                                                          \
        const tmp<CarrierEqn<Type>>& tEqn,                                     \
        const dimensionedScalar& s                                             \
    );

#define CARRIER_SCALAR_EQN_OPERATOR(Op)                                        \
                                                                               \
    template<class Type>                                                       \
    tmp<CarrierEqn<Type>> operator Op                                          \
    (                                                                          \
        const dimensionedScalar& s,                                            \
        const CarrierEqn<Type>& eqn                                            \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    tmp<CarrierEqn<Type>> operator Op                                          \
    (                                                                          \
        const dimensionedScalar& s,                                            \
        const tmp<CarrierEqn<Type>>& tEqn                                      \
    );

//- Multiplication
CARRIER_EQN_SCALAR_OPERATOR(*)
CARRIER_SCALAR_EQN_OPERATOR(*)

//- Division
CARRIER_EQN_SCALAR_OPERATOR(/)

#undef CARRIER_EQN_SCALAR_OPERATOR
#undef CARRIER_SCALAR_EQN_OPERATOR


//- Add to a finite-volume equation
template<class Type>
void operator+=(fvMatrix<Type>& fvEqn, const CarrierEqn<Type>& cEqn);

//- Add to a finite-volume equation
template<class Type>
void operator+=(fvMatrix<Type>& fvEqn, const tmp<CarrierEqn<Type>>& tcEqn);

//- Subtract from a finite-volume equation
template<class Type>
void operator-=(fvMatrix<Type>& fvEqn, const CarrierEqn<Type>& cEqn);

//- Subtract from a finite-volume equation
template<class Type>
void operator-=(fvMatrix<Type>& fvEqn, const tmp<CarrierEqn<Type>>& tcEqn);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CarrierEqn.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
