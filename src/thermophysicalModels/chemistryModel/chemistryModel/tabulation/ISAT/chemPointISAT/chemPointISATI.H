/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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

inline const Foam::chemistryTabulationMethods::ISAT&
Foam::chemPointISAT::table() const
{
    return table_;
}


inline Foam::label Foam::chemPointISAT::nGrowth() const
{
    return nGrowth_;
}


inline Foam::label& Foam::chemPointISAT::completeSpaceSize()
{
    return completeSpaceSize_;
}


inline const Foam::label& Foam::chemPointISAT::completeSpaceSize() const
{
    return completeSpaceSize_;
}


inline const Foam::scalarField& Foam::chemPointISAT::phi() const
{
    return phi_;
}


inline const Foam::scalarField& Foam::chemPointISAT::Rphi() const
{
    return Rphi_;
}


inline const Foam::scalarField& Foam::chemPointISAT::scaleFactor() const
{
    return scaleFactor_;
}


inline const Foam::scalar& Foam::chemPointISAT::tolerance() const
{
    return tolerance_;
}


inline void Foam::chemPointISAT::changeTolerance(scalar newTol)
{
    tolerance_ = newTol;
}


inline Foam::binaryNode*& Foam::chemPointISAT::node()
{
    return node_;
}


inline const Foam::scalarSquareMatrix& Foam::chemPointISAT::A() const
{
    return A_;
}


inline Foam::scalarSquareMatrix& Foam::chemPointISAT::A()
{
    return A_;
}


inline const Foam::scalarSquareMatrix& Foam::chemPointISAT::LT() const
{
    return LT_;
}


inline Foam::scalarSquareMatrix& Foam::chemPointISAT::LT()
{
    return LT_;
}


inline Foam::label Foam::chemPointISAT::nActive() const
{
    return nActive_;
}


inline const Foam::List<Foam::label>&
Foam::chemPointISAT::completeToSimplifiedIndex() const
{
    return completeToSimplifiedIndex_;
}


inline const Foam::List<Foam::label>&
Foam::chemPointISAT::simplifiedToCompleteIndex() const
{
    return simplifiedToCompleteIndex_;
}


inline void Foam::chemPointISAT::increaseNumRetrieve()
{
    this->numRetrieve_++;
}


inline void Foam::chemPointISAT::resetNumRetrieve()
{
    this->numRetrieve_ = 0;
}


inline void Foam::chemPointISAT::increaseNLifeTime()
{
    this->nLifeTime_++;
}


inline Foam::label Foam::chemPointISAT::simplifiedToCompleteIndex
(
    const label i
)
{
    if (i < nActive_)
    {
        return simplifiedToCompleteIndex_[i];
    }
    else if (i == nActive_)
    {
        return completeSpaceSize_ - 3;
    }
    else if (i == nActive_ + 1)
    {
        return completeSpaceSize_ - 2;
    }
    else if (i == nActive_ + 2)
    {
        return completeSpaceSize_ - 1;
    }
    else
    {
        return -1;
    }
}


inline const Foam::label& Foam::chemPointISAT::timeTag()
{
    return timeTag_;
}


inline Foam::label& Foam::chemPointISAT::lastTimeUsed()
{
    return lastTimeUsed_;
}


inline bool& Foam::chemPointISAT::toRemove()
{
    return toRemove_;
}


inline Foam::label& Foam::chemPointISAT::maxNumNewDim()
{
    return maxNumNewDim_;
}


inline const Foam::label& Foam::chemPointISAT::numRetrieve()
{
    return numRetrieve_;
}


inline const Foam::label& Foam::chemPointISAT::nLifeTime()
{
    return nLifeTime_;
}


// ************************************************************************* //
