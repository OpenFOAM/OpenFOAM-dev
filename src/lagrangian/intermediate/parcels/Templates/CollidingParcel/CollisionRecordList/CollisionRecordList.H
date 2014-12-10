/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::CollisionRecordList

Description

SourceFiles
    CollisionRecordListI.H
    CollisionRecordList.C
    CollisionRecordListIO.C

\*---------------------------------------------------------------------------*/

#ifndef CollisionRecordList_H
#define CollisionRecordList_H

#include "DynamicList.H"
#include "PairCollisionRecord.H"
#include "WallCollisionRecord.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class PairType, class WallType>
class CollisionRecordList;

template<class PairType, class WallType>
inline bool operator==
(
    const CollisionRecordList<PairType, WallType>&,
    const CollisionRecordList<PairType, WallType>&
);

template<class PairType, class WallType>
inline bool operator!=
(
    const CollisionRecordList<PairType, WallType>&,
    const CollisionRecordList<PairType, WallType>&
);

template<class PairType, class WallType>
Istream& operator>>(Istream&, CollisionRecordList<PairType, WallType>&);

template<class PairType, class WallType>
Ostream& operator<<(Ostream&, const CollisionRecordList<PairType, WallType>&);


/*---------------------------------------------------------------------------*\
                     Class CollisionRecordList Declaration
\*---------------------------------------------------------------------------*/

template<class PairType, class WallType>
class CollisionRecordList
{
    // Private data

        //- List of active pair collisions
        DynamicList<PairCollisionRecord<PairType> > pairRecords_;

        //- List of active wall collisions
        DynamicList<WallCollisionRecord<WallType> > wallRecords_;


public:

    // Constructors

        //- Construct null
        CollisionRecordList();

        //- Construct from Istream
        CollisionRecordList(Istream&);

        //- Construct from component fields (for IO)
        CollisionRecordList
        (
            const labelField& pairAccessed,
            const labelField& pairOrigProcOfOther,
            const labelField& pairOrigIdOfOther,
            const Field<PairType>& pairData,
            const labelField& wallAccessed,
            const vectorField& wallPRel,
            const Field<WallType>& wallData
        );

    //- Destructor
    ~CollisionRecordList();


    // Member Functions

        //- Return the active pair collisions
        inline const DynamicList<PairCollisionRecord<PairType> >&
            pairRecords() const;

        //- Return the active wall collisions
        inline const DynamicList<WallCollisionRecord<WallType> >&
            wallRecords() const;

        // Fields representing the data from each record, i.e if the
        // records 0-N containing each data members {a, b, c, d...}
        // are organised:
        //
        // a0 b0 c0 d0 ...
        // a1 b1 c1 d1 ...
        // a2 b2 c2 d2 ...
        // ...
        // aN bN cN dN ...
        //
        // Then these field return, for example, (c0, c1, c2,... cN)

        //- Return field of pair accessed from each record, used for
        //  field IO
        labelField pairAccessed() const;

        //- Return field of pair origProcOfOther from each record,
        //  used for field IO
        labelField pairOrigProcOfOther() const;

        //- Return field of pair origIdOfOther from each record, used
        //  for field IO
        labelField pairOrigIdOfOther() const;

        //- Return field of pair data from each record, used for field IO
        Field<PairType> pairData() const;

        //- Return field of wall accessed from each record, used for field IO
        labelField wallAccessed() const;

        //- Return field of wall pRel from each record, used for field IO
        vectorField wallPRel() const;

        //- Return field of wall data from each record, used for field IO
        Field<WallType> wallData() const;

        //- Enquires if the proc and id pair of the other particle are
        //  present in the records.  If so, return non-const access to
        //  the PairCollisionRecord (hence the data) and mark the
        //  PairCollisionRecord as accessed this step, if not, create
        //  the record and return access to it.
        PairCollisionRecord<PairType>& matchPairRecord
        (
            label origProcOfOther,
            label origIdOfOther
        );

        //- Enquire if the specified record exists without modifying
        //  its accessed status
        bool checkPairRecord(label origProcOfOther, label origIdOfOther);

        //- Enquires if the position of wall impact relative to the
        //  particle centre is present in the records.  If so, return
        //  access to the WallCollisionRecord (hence the data) and
        //  mark the WallCollisionRecord as accesses this step, if
        //  not, create the record and return access to it.
        WallCollisionRecord<WallType>& matchWallRecord
        (
            const vector& pRel,
            scalar radius
        );

        //- Enquire if the specified record exists without modifying
        //  its accessed status
        bool checkWallRecord(const vector& pRel, scalar radius);

        //- Update the collision records, deleting any records not
        //  marked as having been accessed, then mark all records as
        //  not accessed ready for the next evaluation
        void update();


        // Member Operators

            void operator=(const CollisionRecordList&);


        // Friend Operators

            friend bool operator== <PairType, WallType>
            (
                const CollisionRecordList<PairType, WallType>&,
                const CollisionRecordList<PairType, WallType>&
            );

            friend bool operator!= <PairType, WallType>
            (
                const CollisionRecordList<PairType, WallType>&,
                const CollisionRecordList<PairType, WallType>&
            );


        // IOstream Operators

            friend Istream& operator>> <PairType, WallType>
            (
                Istream&,
                CollisionRecordList<PairType, WallType>&
            );

            friend Ostream& operator<< <PairType, WallType>
            (
                Ostream&,
                const CollisionRecordList<PairType, WallType>&
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CollisionRecordListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CollisionRecordList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
