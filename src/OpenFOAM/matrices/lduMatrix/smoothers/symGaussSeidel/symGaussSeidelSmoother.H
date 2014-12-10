/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::symGaussSeidelSmoother

Description
    A lduMatrix::smoother for symmetric Gauss-Seidel

SourceFiles
    symGaussSeidelSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef symGaussSeidelSmoother_H
#define symGaussSeidelSmoother_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class symGaussSeidelSmoother Declaration
\*---------------------------------------------------------------------------*/

class symGaussSeidelSmoother
:
    public lduMatrix::smoother
{

public:

    //- Runtime type information
    TypeName("symGaussSeidel");


    // Constructors

        //- Construct from components
        symGaussSeidelSmoother
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );


    // Member Functions

        //- Smooth for the given number of sweeps
        static void smooth
        (
            const word& fieldName,
            scalarField& psi,
            const lduMatrix& matrix,
            const scalarField& source,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const direction cmpt,
            const label nSweeps
        );


        //- Smooth the solution for a given number of sweeps
        virtual void smooth
        (
            scalarField& psi,
            const scalarField& Source,
            const direction cmpt,
            const label nSweeps
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
