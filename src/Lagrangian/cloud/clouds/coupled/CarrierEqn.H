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
{
public:

    // Public Data

        //- Reference to the field
        const VolField<Type>& psi;

        //- Explicit coefficient
        CarrierCoeff<Type, false> Su;

        //- Implicit coefficient
        CarrierCoeff<scalar, true> Sp;


    // Constructors

        //- Construct from a reference to the field
        CarrierEqn(const VolField<Type>& psi);


    // Member Functions

        //- Clear the coefficient values
        void clear();


    // Member Operators

        //- Addition assignment
        void operator+=(const LagrangianEqn<Type>& lEqn);

        //- Subtraction assignment
        void operator-=(const LagrangianEqn<Type>& lEqn);
};


//- Add to a finite-volume equation
template<class Type>
void operator+=(fvMatrix<Type>& fvEqn, const CarrierEqn<Type>& cEqn);

//- Subtract from a finite-volume equation
template<class Type>
void operator-=(fvMatrix<Type>& fvEqn, const CarrierEqn<Type>& cEqn);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CarrierEqn.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
