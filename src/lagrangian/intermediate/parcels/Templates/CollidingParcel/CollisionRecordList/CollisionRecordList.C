/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "CollisionRecordList.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class PairType, class WallType>
Foam::CollisionRecordList<PairType, WallType>::CollisionRecordList()
:
    pairRecords_(),
    wallRecords_()
{}


template<class PairType, class WallType>
Foam::CollisionRecordList<PairType, WallType>::CollisionRecordList(Istream& is)
:
    pairRecords_(is),
    wallRecords_(is)
{
    // Check state of Istream
    is.check
    (
        "Foam::CollisionRecordList<PairType, WallType>::"
        "CollisionRecordList(Foam::Istream&)"
    );
}


template<class PairType, class WallType>
Foam::CollisionRecordList<PairType, WallType>::CollisionRecordList
(
    const labelField& pairAccessed,
    const labelField& pairOrigProcOfOther,
    const labelField& pairOrigIdOfOther,
    const Field<PairType>& pairData,
    const labelField& wallAccessed,
    const vectorField& wallPRel,
    const Field<WallType>& wallData
)
:
    pairRecords_(),
    wallRecords_()
{
    label nPair = pairAccessed.size();

    if
    (
        pairOrigProcOfOther.size() != nPair
     || pairOrigIdOfOther.size() != nPair
     || pairData.size() != nPair
    )
    {
        FatalErrorIn
        (
            "Foam::CollisionRecordList<PairType, WallType>::CollisionRecordList"
            "("
                "const labelField& pairAccessed,"
                "const labelField& pairOrigProcOfOther,"
                "const labelField& pairOrigIdOfOther,"
                "const Field<PairType>& pairData,"
                "const labelField& wallAccessed,"
                "const vectorField& wallPRel,"
                "const Field<WallType>& wallData"
            ")"
        )
            << "Pair field size mismatch." << nl
            << pairAccessed << nl
            << pairOrigProcOfOther << nl
            << pairOrigIdOfOther << nl
            << pairData << nl
            << abort(FatalError);
    }

    forAll(pairAccessed, i)
    {
        pairRecords_.append
        (
            PairCollisionRecord<PairType>
            (
                pairAccessed[i],
                pairOrigProcOfOther[i],
                pairOrigIdOfOther[i],
                pairData[i]
            )
        );
    }

    label nWall = wallAccessed.size();

    if (wallPRel.size() != nWall || wallData.size() != nWall)
    {
        FatalErrorIn
        (
            "Foam::CollisionRecordList<PairType, WallType>::CollisionRecordList"
            "("
                "const labelField& pairAccessed,"
                "const labelField& pairOrigProcOfOther,"
                "const labelField& pairOrigIdOfOther,"
                "const Field<PairType>& pairData,"
                "const labelField& wallAccessed,"
                "const vectorField& wallPRel,"
                "const Field<WallType>& wallData"
            ")"
        )
            << "Wall field size mismatch." << nl
            << wallAccessed << nl
            << wallPRel << nl
            << wallData << nl
            << abort(FatalError);
    }

    forAll(wallAccessed, i)
    {
        wallRecords_.append
        (
            WallCollisionRecord<WallType>
            (
                wallAccessed[i],
                wallPRel[i],
                wallData[i]
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * /

template<class PairType, class WallType>
Foam::CollisionRecordList<PairType, WallType>::~CollisionRecordList()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class PairType, class WallType>
Foam::labelField
Foam::CollisionRecordList<PairType, WallType>::pairAccessed() const
{
    labelField f(pairRecords_.size());

    forAll(pairRecords_, i)
    {
        f[i] = pairRecords_[i].accessed();
    }

    return f;
}


template<class PairType, class WallType>
Foam::labelField
Foam::CollisionRecordList<PairType, WallType>::pairOrigProcOfOther() const
{
    labelField f(pairRecords_.size());

    forAll(pairRecords_, i)
    {
        f[i] = pairRecords_[i].origProcOfOther();
    }

    return f;
}


template<class PairType, class WallType>
Foam::labelField
Foam::CollisionRecordList<PairType, WallType>::pairOrigIdOfOther() const
{
    labelField f(pairRecords_.size());

    forAll(pairRecords_, i)
    {
        f[i] = pairRecords_[i].origIdOfOther();
    }

    return f;
}


template<class PairType, class WallType>
Foam::Field<PairType>
Foam::CollisionRecordList<PairType, WallType>::pairData() const
{
    Field<PairType> f(pairRecords_.size());

    forAll(pairRecords_, i)
    {
        f[i] = pairRecords_[i].collisionData();
    }

    return f;
}


template<class PairType, class WallType>
Foam::labelField
Foam::CollisionRecordList<PairType, WallType>::wallAccessed() const
{
    labelField f(wallRecords_.size());

    forAll(wallRecords_, i)
    {
        f[i] = wallRecords_[i].accessed();
    }

    return f;
}


template<class PairType, class WallType>
Foam::vectorField
Foam::CollisionRecordList<PairType, WallType>::wallPRel() const
{
    vectorField f(wallRecords_.size());

    forAll(wallRecords_, i)
    {
        f[i] = wallRecords_[i].pRel();
    }

    return f;
}


template<class PairType, class WallType>
Foam::Field<WallType>
Foam::CollisionRecordList<PairType, WallType>::wallData() const
{
    Field<WallType> f(wallRecords_.size());

    forAll(wallRecords_, i)
    {
        f[i] = wallRecords_[i].collisionData();
    }

    return f;
}


template<class PairType, class WallType>
Foam::PairCollisionRecord<PairType>&
Foam::CollisionRecordList<PairType, WallType>::matchPairRecord
(
    label origProcOfOther,
    label origIdOfOther
)
{
    // Returning the first record that matches the particle
    // identifiers.  Two records with the same identification is not
    // supported.

    forAll(pairRecords_, i)
    {
        PairCollisionRecord<PairType>& pCR = pairRecords_[i];

        if (pCR.match(origProcOfOther, origIdOfOther))
        {
            pCR.setAccessed();

            return pCR;
        }
    }

    // Record not found, create a new one and return it as the last
    // member of the list.  Setting the status of the record to be accessed
    // on construction.

    pairRecords_.append
    (
        PairCollisionRecord<PairType>(true, origProcOfOther, origIdOfOther)
    );

    return pairRecords_.last();
}


template<class PairType, class WallType>
bool Foam::CollisionRecordList<PairType, WallType>::checkPairRecord
(
    label origProcOfOther,
    label origIdOfOther
)
{
    forAll(pairRecords_, i)
    {
        PairCollisionRecord<PairType>& pCR = pairRecords_[i];

        if (pCR.match(origProcOfOther, origIdOfOther))
        {
            return true;
        }
    }

    return false;
}


template<class PairType, class WallType>
Foam::WallCollisionRecord<WallType>&
Foam::CollisionRecordList<PairType, WallType>::matchWallRecord
(
    const vector& pRel,
    scalar radius
)
{
    // Returning the first record that matches the relative position.
    // Two records with the same relative position is not supported.

    forAll(wallRecords_, i)
    {
        WallCollisionRecord<WallType>& wCR = wallRecords_[i];

        if (wCR.match(pRel, radius))
        {
            wCR.setAccessed();

            return wCR;
        }
    }

    // Record not found, create a new one and return it as the last
    // member of the list.  Setting the status of the record to be accessed
    // on construction.

    wallRecords_.append(WallCollisionRecord<WallType>(true, pRel));

    return wallRecords_.last();
}


template<class PairType, class WallType>
bool Foam::CollisionRecordList<PairType, WallType>::checkWallRecord
(
    const vector& pRel,
    scalar radius
)
{
    forAll(wallRecords_, i)
    {
        WallCollisionRecord<WallType>& wCR = wallRecords_[i];

        if (wCR.match(pRel, radius))
        {
            return true;
        }
    }

    return false;
}


template<class PairType, class WallType>
void Foam::CollisionRecordList<PairType, WallType>::update()
{
    {
        DynamicList<PairCollisionRecord<PairType> > updatedRecords;

        forAll(pairRecords_, i)
        {
            if (pairRecords_[i].accessed())
            {
                pairRecords_[i].setUnaccessed();

                updatedRecords.append(pairRecords_[i]);
            }
        }

        pairRecords_ = updatedRecords;
    }

    {
        DynamicList<WallCollisionRecord<WallType> > updatedRecords;

        forAll(wallRecords_, i)
        {
            if (wallRecords_[i].accessed())
            {
                wallRecords_[i].setUnaccessed();

                updatedRecords.append(wallRecords_[i]);
            }
        }

        wallRecords_ = updatedRecords;
    }
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class PairType, class WallType>
void Foam::CollisionRecordList<PairType, WallType>::operator=
(
    const CollisionRecordList<PairType, WallType>& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "Foam::CollisionRecordList<PairType, WallType>::operator="
            "(const Foam::CollisionRecordList<PairType, WallType>&)"
        )
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    pairRecords_ = rhs.pairRecords_;
    wallRecords_ = rhs.wallRecords_;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

template<class PairType, class WallType>
inline bool Foam::operator==
(
    const CollisionRecordList<PairType, WallType>& a,
    const CollisionRecordList<PairType, WallType>& b
)
{
    return
    (
        a.pairRecords_ == b.pairRecords_
     && a.wallRecords_ == b.wallRecords_
    );
}


template<class PairType, class WallType>
inline bool Foam::operator!=
(
    const CollisionRecordList<PairType, WallType>& a,
    const CollisionRecordList<PairType, WallType>& b
)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class PairType, class WallType>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    CollisionRecordList<PairType, WallType>& cRL
)
{
    is  >> cRL.pairRecords_ >> cRL.wallRecords_;

    // Check state of Istream
    is.check
    (
        "Foam::Istream& Foam::operator>>"
        "(Foam::Istream&, Foam::CollisionRecordList<PairType, WallType>&)"
    );

    return is;
}


template<class PairType, class WallType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const CollisionRecordList<PairType, WallType>& cRL
)
{
    os  << cRL.pairRecords_ << cRL.wallRecords_;

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<(Foam::Ostream&, "
        "const Foam::CollisionRecordList<PairType, WallType>&)"
    );

    return os;
}


// ************************************************************************* //
