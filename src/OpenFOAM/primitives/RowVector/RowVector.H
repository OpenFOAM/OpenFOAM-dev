/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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
    Foam::RowVector

Description
    Templated 3D row-vector derived from MatrixSpace adding construction from
    3 components and element access using x(), y() and z().

SourceFiles
    RowVectorI.H

\*---------------------------------------------------------------------------*/

#ifndef RowVector_H
#define RowVector_H

#include "MatrixSpace.H"
#include "Vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class RowVector Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class RowVector
:
    public MatrixSpace<RowVector<Cmpt>, Cmpt, 1, 3>
{

public:

    //- Equivalent type of labels used for valid component indexing
    typedef RowVector<label> labelType;


    //- Component labeling enumeration
    enum components { X, Y, Z };


    // Constructors

        //- Construct null
        inline RowVector();

        //- Construct initialised to zero
        inline RowVector(const Foam::zero);

        //- Construct given VectorSpace of the same rank
        template<class Cmpt2>
        inline RowVector(const MatrixSpace<RowVector<Cmpt2>, Cmpt2, 1, 3>&);

        //- Construct given three components
        inline RowVector(const Cmpt& rvx, const Cmpt& rvy, const Cmpt& rvz);

        //- Construct from Istream
        inline RowVector(Istream&);


    // Member Functions

        // Access

            inline const Cmpt& x() const;
            inline const Cmpt& y() const;
            inline const Cmpt& z() const;

            inline Cmpt& x();
            inline Cmpt& y();
            inline Cmpt& z();
};


template<class Cmpt>
class typeOfTranspose<Cmpt, Vector<Cmpt>>
{
public:

    typedef RowVector<Cmpt> type;
};


template<class Cmpt>
class typeOfTranspose<Cmpt, RowVector<Cmpt>>
{
public:

    typedef Vector<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "RowVectorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
