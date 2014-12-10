/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::LUscalarMatrix

Description
    Foam::LUscalarMatrix

SourceFiles
    LUscalarMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef LUscalarMatrix_H
#define LUscalarMatrix_H

#include "scalarMatrices.H"
#include "labelList.H"
#include "FieldField.H"
#include "lduInterfaceFieldPtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lduMatrix;
class procLduMatrix;

/*---------------------------------------------------------------------------*\
                           Class LUscalarMatrix Declaration
\*---------------------------------------------------------------------------*/

class LUscalarMatrix
:
    public scalarSquareMatrix
{
    // Private data

        //- Communicator to use
        const label comm_;

        //- Processor matrix offsets
        labelList procOffsets_;

        //- The pivot indices used in the LU decomposition
        labelList pivotIndices_;

        //- Convert the given lduMatrix into this LUscalarMatrix
        void convert
        (
            const lduMatrix& ldum,
            const FieldField<Field, scalar>& interfaceCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );

        //- Convert the given list of procLduMatrix into this LUscalarMatrix
        //  on the master processor
        void convert(const PtrList<procLduMatrix>& lduMatrices);


        //- Print the ratio of the mag-sum of the off-diagonal coefficients
        //  to the mag-diagonal
        void printDiagonalDominance() const;


public:

    // Declare name of the class and its debug switch
    ClassName("LUscalarMatrix");

    // Constructors

        //- Construct from scalarSquareMatrix and perform LU decomposition
        LUscalarMatrix(const scalarSquareMatrix&);

        //- Construct from lduMatrix and perform LU decomposition
        LUscalarMatrix
        (
            const lduMatrix&,
            const FieldField<Field, scalar>& interfaceCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );


    // Member Functions

        //- Solve the matrix using the LU decomposition with pivoting
        //  returning the solution in the source
        template<class T>
        void solve(Field<T>& source) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "LUscalarMatrixTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
