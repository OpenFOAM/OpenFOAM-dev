/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "cutPolyValue.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Type Foam::cutPoly::edgeCutValue
(
    const edge& e,
    const scalar lambda,
    const Field<Type>& pPsis
)
{
    return (1 - lambda)*pPsis[e[0]] + lambda*pPsis[e[1]];
}


template<class Type>
Type Foam::cutPoly::edgeCutValue
(
    const edge& e,
    const scalarField& pAlphas,
    const scalar isoAlpha,
    const Field<Type>& pPsis
)
{
    return edgeCutValue(e, edgeCutLambda(e, pAlphas, isoAlpha), pPsis);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
const Foam::Pair<Type>&
Foam::cutPoly::FaceCutValues<Type>::const_iterator::cutPsis(const label i) const
{
    if (i != iAndCutPsis_.first())
    {
        const face& f = fValues_.f_;
        const labelPair& fCut = fValues_.fCuts_[i];

        iAndCutPsis_.first() = i;

        forAll(iAndCutPsis_.second(), j)
        {
            iAndCutPsis_.second()[j] =
                edgeCutValue
                (
                    f.faceEdge(fCut[j]),
                    fValues_.pAlphas_,
                    fValues_.isoAlpha_,
                    fValues_.pPsis_
                );
        }
    }

    return iAndCutPsis_.second();
}


template<class Type>
Foam::label
Foam::cutPoly::FaceCutValues<Type>::const_iterator::size(const label i) const
{
    const face& f = fValues_.f_;
    const List<labelPair>& fCuts = fValues_.fCuts_;
    const labelPair& fCut0 = fCuts[i];
    const labelPair& fCut1 = fCuts[(i + 1) % fCuts.size()];

    const label dfCut =
        fValues_.below_ == separatedBelow
      ? fCut0[1] - fCut0[0]
      : fCut1[0] - fCut0[1];

    return 2 + ((dfCut + f.size()) % f.size());
}


template<class Type>
const Type Foam::cutPoly::FaceCutValues<Type>::const_iterator::psi
(
    const label i,
    const label j
) const
{
    const face& f = fValues_.f_;
    const labelPair& fCut = fValues_.fCuts_[i];

    if (j < 2)
    {
        const label fCutj = fValues_.below_ == separatedBelow ? 1 - j : j;
        return cutPsis(i)[fCutj];
    }
    else
    {
        const label fCutj = fValues_.below_ == separatedBelow ? 0 : 1;
        return fValues_.pPsis_[f[(j - 2 + fCut[fCutj] + 1) % f.size()]];
    }
}


template<class Type>
Type Foam::cutPoly::CellCutValues<Type>::const_iterator::psi
(
    const label i,
    const label j
)
{
    const label iStar = i % cValues_.cCuts_.size();
    const label jStar = j % cValues_.cCuts_[iStar].size();

    const label cei = cValues_.cCuts_[iStar][jStar];
    const label cfi = cValues_.cAddr_.ceiToCfiAndFei()[cei][0][0];
    const label fei = cValues_.cAddr_.ceiToCfiAndFei()[cei][0][1];

    return
        edgeCutValue
        (
            cValues_.fs_[cValues_.c_[cfi]].faceEdge(fei),
            cValues_.pAlphas_,
            cValues_.isoAlpha_,
            cValues_.pPsis_
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cutPoly::FaceCutValues<Type>::const_iterator::const_iterator
(
    const FaceCutValues<Type>& fValues,
    const label i,
    const label j
)
:
    fValues_(fValues),
    i_(i),
    j_(j),
    iAndCutPsis_(-1, Pair<Type>(Zero, Zero))
{}


template<class Type>
Foam::cutPoly::FaceCutValues<Type>::FaceCutValues
(
    const face& f,
    const List<labelPair>& fCuts,
    const Field<Type>& pPsis,
    const scalarField& pAlphas,
    const scalar isoAlpha,
    const bool below
)
:
    f_(f),
    fCuts_(fCuts),
    pPsis_(pPsis),
    pAlphas_(pAlphas),
    isoAlpha_(isoAlpha),
    below_(below)
{}


template<class Type>
Foam::cutPoly::CellCutValues<Type>::const_iterator::const_iterator
(
    const CellCutValues<Type>& cValues,
    const label i,
    const label j
)
:
    cValues_(cValues),
    i_(i),
    j_(j),
    psis_(psi(i, j), psi(i, j + 1))
{}


template<class Type>
Foam::cutPoly::CellCutValues<Type>::CellCutValues
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const labelListList& cCuts,
    const faceList& fs,
    const Field<Type>& pPsis,
    const scalarField& pAlphas,
    const scalar isoAlpha
)
:
    c_(c),
    cAddr_(cAddr),
    cCuts_(cCuts),
    fs_(fs),
    pPsis_(pPsis),
    pAlphas_(pAlphas),
    isoAlpha_(isoAlpha)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Type Foam::cutPoly::FaceCutValues<Type>::const_iterator::next() const
{
    const label j = (j_ + 1) % size(i_);
    const label i = i_ + (fValues_.below_ != separatedBelow && j == 0);
    return i == fValues_.fCuts_.size() ? psi(0, 0) : psi(i, j);
}


template<class Type>
typename Foam::cutPoly::FaceCutValues<Type>::const_iterator
Foam::cutPoly::FaceCutValues<Type>::begin() const
{
    return const_iterator(*this, 0, 0);
}


template<class Type>
typename Foam::cutPoly::FaceCutValues<Type>::const_iterator
Foam::cutPoly::FaceCutValues<Type>::end() const
{
    return const_iterator(*this, this->fCuts_.size(), 0);
}


template<class Type>
const Type& Foam::cutPoly::CellCutValues<Type>::const_iterator::next() const
{
    return psis_.second();
}


template<class Type>
typename Foam::cutPoly::CellCutValues<Type>::const_iterator
Foam::cutPoly::CellCutValues<Type>::begin() const
{
    return const_iterator(*this, 0, 0);
}


template<class Type>
typename Foam::cutPoly::CellCutValues<Type>::const_iterator
Foam::cutPoly::CellCutValues<Type>::end() const
{
    return const_iterator(*this, this->cCuts_.size(), 0);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
bool Foam::cutPoly::FaceCutValues<Type>::const_iterator::operator==
(
    const const_iterator& it
) const
{
    return it.i_ == i_ && it.j_ == j_;
}


template<class Type>
bool Foam::cutPoly::FaceCutValues<Type>::const_iterator::operator!=
(
    const const_iterator& it
) const
{
    return !(it == *this);
}


template<class Type>
Type Foam::cutPoly::FaceCutValues<Type>::const_iterator::operator*() const
{
    return psi(i_, j_);
}


template<class Type>
inline typename Foam::cutPoly::FaceCutValues<Type>::const_iterator&
Foam::cutPoly::FaceCutValues<Type>::const_iterator::operator++()
{
    j_ = (j_ + 1) % size(i_);
    i_ += j_ == 0;
    return *this;
}


template<class Type>
inline typename Foam::cutPoly::FaceCutValues<Type>::const_iterator
Foam::cutPoly::FaceCutValues<Type>::const_iterator::operator++(int)
{
    const const_iterator it(*this);
    ++ *this;
    return it;
}


template<class Type>
bool Foam::cutPoly::CellCutValues<Type>::const_iterator::operator==
(
    const const_iterator& it
) const
{
    return it.i_ == i_ && it.j_ == j_;
}

template<class Type>
bool Foam::cutPoly::CellCutValues<Type>::const_iterator::operator!=
(
    const const_iterator& it
) const
{
    return !(it == *this);
}


template<class Type>
const Type&
Foam::cutPoly::CellCutValues<Type>::const_iterator::operator*() const
{
    return psis_.first();
}


template<class Type>
inline typename Foam::cutPoly::CellCutValues<Type>::const_iterator&
Foam::cutPoly::CellCutValues<Type>::const_iterator::operator++()
{
    ++ j_;

    if (j_ == cValues_.cCuts_[i_].size())
    {
        j_ = 0;
        ++ i_;
        psis_.first() = psi(i_, j_);
    }
    else
    {
        psis_.first() = psis_.second();
    }

    psis_.second() = psi(i_, j_ + 1);

    return *this;
}


template<class Type>
inline typename Foam::cutPoly::CellCutValues<Type>::const_iterator
Foam::cutPoly::CellCutValues<Type>::const_iterator::operator++(int)
{
    const const_iterator it(*this);
    ++ *this;
    return it;
}


// ************************************************************************* //
