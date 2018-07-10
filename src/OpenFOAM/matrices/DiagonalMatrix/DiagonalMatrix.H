/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::DiagonalMatrix

Description
    DiagonalMatrix<Type> is a 2D diagonal matrix of objects
    of type Type, size nxn

SourceFiles
    DiagonalMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef DiagonalMatrix_H
#define DiagonalMatrix_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * *  * * * * * * Class Forward declaration  * * * * * * * * * * * //

template<class Form, class Type> class Matrix;

/*---------------------------------------------------------------------------*\
                           Class DiagonalMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class DiagonalMatrix
:
    public List<Type>
{
public:

    // Constructors

        //- Null constructor.
        DiagonalMatrix<Type>();

        //- Construct from diagonal component of a Matrix
        template<class Form>
        DiagonalMatrix<Type>(const Matrix<Form, Type>&);

        //- Construct empty from size
        DiagonalMatrix<Type>(const label size);

        //- Construct from size and a value
        DiagonalMatrix<Type>(const label, const Type&);


    // Member functions

        //- Invert the diagonal matrix and return itself
        DiagonalMatrix<Type>& invert();
};


// Global functions

//- Return the diagonal Matrix inverse
template<class Type>
DiagonalMatrix<Type> inv(const DiagonalMatrix<Type>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DiagonalMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
