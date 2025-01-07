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
    Foam::CarrierCoeff

Description
    Class to store a coefficient of a carrier equation

SourceFiles
    CarrierCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef CarrierCoeff_H
#define CarrierCoeff_H

#include "DimensionedField.H"
#include "volMesh.H"
#include "LagrangianCoeff.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class CarrierCoeff Declaration
\*---------------------------------------------------------------------------*/

template<class Type, bool Implicit>
class CarrierCoeff
{
    // Private Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Coefficient
        autoPtr<DimensionedField<Type, volMesh>> S_;


    // Private Member Functions

        //- Construct the coefficient field for a given Lagrangian field
        void initialise(const LagrangianSubField<Type>& lField);


public:

    // Constructors

        //- Construct from a reference to the field
        template<class PsiType>
        CarrierCoeff(const VolField<PsiType>& psi);

        //- Construct as copy
        CarrierCoeff(const CarrierCoeff<Type, Implicit>& coeff);

        //- Construct as copy or re-use as specified
        CarrierCoeff(CarrierCoeff<Type, Implicit>& coeff, const bool reuse);

        //- Move construct
        CarrierCoeff(CarrierCoeff<Type, Implicit>&& coeff);


    // Member Functions

        //- Determine whether this coefficient has values or not
        bool valid() const;

        //- Clear the coefficient values
        void clear();

        //- Access the field
        const DimensionedField<Type, volMesh>& S() const;

        //- In-place negation
        void negate();


    // Member Operators

        //- Addition assignment
        void operator+=(const LagrangianCoeff<Type, Implicit>& lCoeff);

        //- Subtraction assignment
        void operator-=(const LagrangianCoeff<Type, Implicit>& lCoeff);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CarrierCoeff.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
