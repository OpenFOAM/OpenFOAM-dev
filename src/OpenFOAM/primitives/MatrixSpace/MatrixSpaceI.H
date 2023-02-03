/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include <type_traits>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::MatrixSpace()
{}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::MatrixSpace
(
    const Foam::zero
)
:
    MatrixSpace::vsType(Zero)
{}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
template<class Form2, class Cmpt2>
inline Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::MatrixSpace
(
    const VectorSpace<Form2, Cmpt2, Mrows*Ncols>& vs
)
:
    MatrixSpace::vsType(vs)
{}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
template
<
    template<class, Foam::direction, Foam::direction> class Block2,
    Foam::direction BRowStart,
    Foam::direction BColStart
>
inline Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::MatrixSpace
(
    const Block2<Form, BRowStart, BColStart>& block
)
{
    for (direction i=0; i<Mrows; ++i)
    {
        for (direction j=0; j<Ncols; ++j)
        {
            operator()(i, j) = block(i, j);
        }
    }
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::MatrixSpace(Istream& is)
:
    MatrixSpace::vsType(is)
{}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
template<class SubTensor, Foam::direction BRowStart, Foam::direction BColStart>
inline Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::
ConstBlock<SubTensor, BRowStart, BColStart>::
ConstBlock(const msType& matrix)
:
    matrix_(matrix)
{
    static_assert
    (
        msType::mRows >= BRowStart + mRows,
        "Rows in block > rows in matrix"
    );
    static_assert
    (
        msType::nCols >= BColStart + nCols,
        "Columns in block > columns in matrix"
    );
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
template<class SubTensor, Foam::direction BRowStart, Foam::direction BColStart>
inline Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::
Block<SubTensor, BRowStart, BColStart>::
Block(msType& matrix)
:
    matrix_(matrix)
{
    static_assert
    (
        msType::mRows >= BRowStart + mRows,
        "Rows in block > rows in matrix"
    );
    static_assert
    (
        msType::nCols >= BColStart + nCols,
        "Columns in block > columns in matrix"
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
template<Foam::direction Row, Foam::direction Col>
inline const Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::elmt() const
{
    static_assert(Row < Mrows && Col < Ncols, "Address outside matrix");
    return this->v_[Row*Ncols + Col];
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
template<Foam::direction Row, Foam::direction Col>
inline Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::elmt()
{
    static_assert(Row < Mrows && Col < Ncols, "Address outside matrix");
    return this->v_[Row*Ncols + Col];
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline const Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::xx() const
{
    return elmt<0, 0>();
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::xx()
{
    return elmt<0, 0>();
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline const Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::xy() const
{
    return elmt<0,1>();
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::xy()
{
    return elmt<0,1>();
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline const Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::xz() const
{
    return elmt<0,2>();
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::xz()
{
    return elmt<0,2>();
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline const Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::yx() const
{
    return elmt<1,0>();
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::yx()
{
    return elmt<1,0>();
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline const Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::yy() const
{
    return elmt<1,1>();
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::yy()
{
    return elmt<1,1>();
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline const Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::yz() const
{
    return elmt<1,2>();
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::yz()
{
    return elmt<1,2>();
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline const Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::zx() const
{
    return elmt<2,0>();
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::zx()
{
    return elmt<2,0>();
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline const Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::zy() const
{
    return elmt<2,1>();
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::zy()
{
    return elmt<2,1>();
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline const Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::zz() const
{
    return elmt<2,2>();
}

template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::zz()
{
    return elmt<2,2>();
}

template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>
Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::identityMap()
{
    static_assert(Mrows == Ncols, "Matrix is not square");
    msType result(Zero);

    for (direction i=0; i<Ncols; ++i)
    {
        result(i, i) = 1;
    }

    return result;
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline typename Foam::typeOfTranspose<Cmpt, Form>::type
Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::T() const
{
    typename typeOfTranspose<Cmpt, Form>::type result;

    for (direction i=0; i<Mrows; ++i)
    {
        for (direction j=0; j<Ncols; ++j)
        {
            result(j,i) = operator()(i, j);
        }
    }

    return result;
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
template
<
    class SubTensor,
    Foam::direction BRowStart,
    Foam::direction BColStart
>
inline typename Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::template
    ConstBlock<SubTensor, BRowStart, BColStart>
Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::block() const
{
    return *this;
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
template
<
    class SubTensor,
    Foam::direction BRowStart,
    Foam::direction BColStart
>
inline
typename Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::template
    Block<SubTensor, BRowStart, BColStart>
Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::block()
{
    return *this;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline const Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::operator()
(
    const direction& i,
    const direction& j
) const
{
    #ifdef FULLDEBUG
    if (i > Mrows-1 || j > Ncols-1)
    {
        FatalErrorInFunction
            << "indices out of range"
            << abort(FatalError);
    }
    #endif

    return this->v_[i*Ncols + j];
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline Cmpt& Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::operator()
(
    const direction& i,
    const direction& j
)
{
    #ifdef FULLDEBUG
    if (i > Mrows-1 || j > Ncols-1)
    {
        FatalErrorInFunction
            << "indices out of range"
            << abort(FatalError);
    }
    #endif

    return this->v_[i*Ncols + j];
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
template<class SubTensor, Foam::direction BRowStart, Foam::direction BColStart>
inline SubTensor
Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::
ConstBlock<SubTensor, BRowStart, BColStart>::
operator()() const
{
    return *this;
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
template<class SubTensor, Foam::direction BRowStart, Foam::direction BColStart>
inline const Cmpt&
Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::
ConstBlock<SubTensor, BRowStart, BColStart>::
operator()(const direction i, const direction j) const
{
    return matrix_(BRowStart + i, BColStart + j);
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
template<class SubTensor, Foam::direction BRowStart, Foam::direction BColStart>
inline SubTensor
Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::
Block<SubTensor, BRowStart, BColStart>::
operator()() const
{
    SubTensor st;

    for (direction i=0; i<SubTensor::mRows; ++i)
    {
        for (direction j=0; j<SubTensor::nCols; ++j)
        {
            st[i*SubTensor::nCols + j] = operator()(i, j);
        }
    }

    return st;
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
template<class SubTensor, Foam::direction BRowStart, Foam::direction BColStart>
inline const Cmpt&
Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::
Block<SubTensor, BRowStart, BColStart>::
operator()(const direction i, const direction j) const
{
    return matrix_(BRowStart + i, BColStart + j);
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
template<class SubTensor, Foam::direction BRowStart, Foam::direction BColStart>
inline Cmpt&
Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::
Block<SubTensor, BRowStart, BColStart>::
operator()(const direction i, const direction j)
{
    return matrix_(BRowStart + i, BColStart + j);
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
inline void Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::operator=
(
    const Foam::zero
)
{
    MatrixSpace::vsType::operator=(Zero);
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
template<class Form2>
inline void Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::operator&=
(
    const MatrixSpace<Form, Cmpt, Ncols, Ncols>& matrix
)
{
    *this = *this & matrix;
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
template
<
    template<class, Foam::direction, Foam::direction> class Block2,
    Foam::direction BRowStart,
    Foam::direction BColStart
>
inline void Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::operator=
(
    const Block2<Form, BRowStart, BColStart>& block
)
{
    for (direction i = 0; i < Mrows; ++i)
    {
        for (direction j = 0; j < Ncols; ++j)
        {
            operator()(i, j) = block(i, j);
        }
    }
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
template<class SubTensor, Foam::direction BRowStart, Foam::direction BColStart>
template<class Form2>
inline void
Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::
Block<SubTensor, BRowStart, BColStart>::
operator=
(
    const MatrixSpace<Form2, Cmpt, SubTensor::mRows, SubTensor::nCols>& matrix
)
{
    for (direction i=0; i<mRows; ++i)
    {
        for (direction j=0; j<nCols; ++j)
        {
            operator()(i,j) = matrix(i,j);
        }
    }
}


template<class Form, class Cmpt, Foam::direction Mrows, Foam::direction Ncols>
template<class SubTensor, Foam::direction BRowStart, Foam::direction BColStart>
template<class VSForm>
inline void
Foam::MatrixSpace<Form, Cmpt, Mrows, Ncols>::
Block<SubTensor, BRowStart, BColStart>::
operator=
(
    const VectorSpace<VSForm, Cmpt, SubTensor::mRows>& v
)
{
    static_assert(nCols == 1, "Matrix must have a single column");

    for (direction i=0; i<SubTensor::mRows; ++i)
    {
        operator()(i,0) = v[i];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

template<class Form, class Cmpt, direction Mrows, direction Ncols>
inline typename typeOfTranspose<Cmpt, Form>::type T
(
    const MatrixSpace<Form, Cmpt, Ncols, Mrows>& matrix
)
{
    return matrix.T();
}


template<class Form, class Cmpt, direction Ncmpts>
inline typename typeOfTranspose<Cmpt, Form>::type T
(
    const VectorSpace<Form, Cmpt, Ncmpts>& v
)
{
    typename typeOfTranspose<Cmpt, Form>::type result;

    for (direction i=0; i<Ncmpts; ++i)
    {
        result[i] = v[i];
    }

    return result;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template
<
    class Form1,
    class Form2,
    class Cmpt,
    direction Mrows1,
    direction Ncols1,
    direction Mrows2,
    direction Ncols2
>
inline typename typeOfInnerProduct<Cmpt, Form1, Form2>::type operator&
(
    const MatrixSpace<Form1, Cmpt, Mrows1, Ncols1>& matrix1,
    const MatrixSpace<Form2, Cmpt, Mrows2, Ncols2>& matrix2
)
{
    static_assert
    (
        Ncols1 == Mrows2,
        "Number of columns in matrix 1 != number of rows in matrix 2"
    );

    typename typeOfInnerProduct<Cmpt, Form1, Form2>::type result(Zero);

    for (direction i=0; i<Mrows1; ++i)
    {
        for (direction j=0; j<Ncols2; ++j)
        {
            for (direction k=0; k<Mrows2; k++)
            {
                result(i, j) += matrix1(i, k)*matrix2(k, j);
            }
        }
    }

    return result;
}


template<class Form, class VSForm, class Cmpt, direction Mrows, direction Ncols>
inline typename typeOfInnerProduct<Cmpt, Form, VSForm>::type operator&
(
    const MatrixSpace<Form, Cmpt, Mrows, Ncols>& matrix,
    const VectorSpace<VSForm, Cmpt, Ncols>& v
)
{
    typename typeOfInnerProduct<Cmpt, Form, VSForm>::type result(Zero);

    for (direction i=0; i<Mrows; ++i)
    {
        for (direction j=0; j<Ncols; ++j)
        {
            result[i] += matrix(i, j)*v[j];
        }
    }

    return result;
}


template
<
    class Form1,
    class Form2,
    class Cmpt,
    direction Ncmpts1,
    direction Ncmpts2
>
inline typename typeOfOuterProduct<Cmpt, Form1, Form2>::type operator*
(
    const VectorSpace<Form1, Cmpt, Ncmpts1>& v1,
    const VectorSpace<Form2, Cmpt, Ncmpts2>& v2
)
{
    typename typeOfOuterProduct<Cmpt, Form1, Form2>::type result;

    for (direction i=0; i<Ncmpts1; ++i)
    {
        for (direction j=0; j<Ncmpts2; ++j)
        {
            result(i, j) = v1[i]*v2[j];
        }
    }

    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
