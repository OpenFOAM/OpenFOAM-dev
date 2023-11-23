/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::regionSolvers::size() const
{
    return solvers_.size();
}


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

inline Foam::regionSolvers::iterator::iterator(regionSolvers& rs)
:
   regionSolvers_(rs),
   index_(0)
{}


inline bool Foam::regionSolvers::iterator::operator==
(
    const iterator& iter
) const
{
    return index_ == iter.index_;
}


inline bool Foam::regionSolvers::iterator::operator!=
(
    const iterator& iter
) const
{
    return !operator==(iter);
}


inline Foam::solver& Foam::regionSolvers::iterator::operator*()
{
    return regionSolvers_.solvers_[index_];
}


inline Foam::solver& Foam::regionSolvers::iterator::operator()()
{
    return regionSolvers_.solvers_[index_];
}


inline Foam::solver* Foam::regionSolvers::iterator::operator->()
{
    return regionSolvers_.solvers_(index_);
}


inline Foam::regionSolvers::iterator
Foam::regionSolvers::iterator::operator++()
{
    if (++index_ >= regionSolvers_.solvers_.size())
    {
        // Set index to -1 to indicate end
        index_ = -1;

        // Reset the prefix to global space padding
        regionSolvers_.setGlobalPrefix();
    }
    else
    {
        // Set the prefix for region corresponding to index_
        regionSolvers_.setPrefix(index_);
    }

    return *this;
}


inline Foam::regionSolvers::iterator
Foam::regionSolvers::iterator::operator++(int)
{
    iterator tmp(*this);
    operator++();
    return tmp;
}


inline Foam::regionSolvers::iterator Foam::regionSolvers::begin()
{
    // Set the prefix for region 0
    setPrefix(0);

    // Return the iterator for region 0
    return iterator(*this);
}


inline Foam::regionSolvers::iterator Foam::regionSolvers::end()
{
    iterator endIter(*this);

    // Set index to -1 to indicate end
    endIter.index_ = -1;

    return endIter;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::regionSolvers::operator PtrList<solver>&()
{
    return solvers_;
}


// ************************************************************************* //
