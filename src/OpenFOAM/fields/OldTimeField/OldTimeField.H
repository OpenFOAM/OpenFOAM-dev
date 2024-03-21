/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::OldTimeField

Description
    Class to add into field types to provide old-time storage and retrieval

SourceFiles
    OldTimeFieldI.H
    OldTimeField.C

\*---------------------------------------------------------------------------*/

#ifndef OldTimeField_H
#define OldTimeField_H

#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class FieldType>
class OldTimeField;


/*---------------------------------------------------------------------------*\
                        Using OldTimeVoid Declaration
\*---------------------------------------------------------------------------*/

template<typename T>
using OldTimeVoid = void;


/*---------------------------------------------------------------------------*\
                    Struct OldTimeBaseFieldType Declaration
\*---------------------------------------------------------------------------*/

template<class FieldType, typename = void>
struct OldTimeBaseFieldType
{
    nil operator()(const OldTimeField<FieldType>& otf) const
    {
        return nil();
    }
};

template<class FieldType>
struct OldTimeBaseFieldType
<
    FieldType,
    OldTimeVoid<typename FieldType::Base::OldTime>
>
{
    typedef typename FieldType::Base::OldTime type;

    const type& operator()(const OldTimeField<FieldType>& otf) const
    {
        return static_cast<const type&>(otf.field());
    }
};


/*---------------------------------------------------------------------------*\
                          Class OldTimeField Declaration
\*---------------------------------------------------------------------------*/

template<class FieldType>
class OldTimeField
{
public:

    // Public Typedefs

        //- The old-time field type
        typedef OldTimeField<FieldType> OldTime;


private:

    // Private Data

        //- Current time index
        mutable label timeIndex_;

        //- The old time field pointer or reference
        mutable tmp<FieldType> tfield0_;


    // Private Member Functions

        //- Get a reference to the field
        const FieldType& field() const;

        //- Get a non-const reference to the field
        FieldType& fieldRef();

        //- Store the old-time fields. Inner recursion.
        void storeOldTimesInner() const;

        //- Set the oldest field pointer to nullObjectPtr. Inner recursion.
        void nullOldestTimeInner();

        //- Set the field reference in the base class. Overload for when the
        //  base class is not an old-time field.
        void setBase(const nil&) const;

        //- Set the field reference in the base class. Overload for when the
        //  base class is an old-time field.
        template<class OldTimeBaseField>
        void setBase(const OldTimeBaseField& otbf) const;

        //- Set the field reference in the base class
        void setBase() const;


protected:

    // Protected Member Functions

        //- Read old-time field from file if it is present
        bool readOldTimeIfPresent();

        //- Copy the old times from the given field
        void copyOldTimes
        (
            const IOobject& io,
            const OldTimeField<FieldType>& otf
        );

        //- Copy the old times from the given field
        void copyOldTimes
        (
            const word& newName,
            const OldTimeField<FieldType>& otf
        );


public:

    //- Declare friendship with other old-time fields
    template<class OtherFieldType>
    friend class OldTimeField;

    //- Declare friendship with the base old-time field casting struct
    template<class OtherFieldType, typename>
    friend struct OldTimeBaseFieldType;


    // Constructors

        //- Construct from a time index
        OldTimeField(const label timeIndex);

        //- Copy construct
        OldTimeField(const OldTimeField<FieldType>& otf);

        //- Move construct
        OldTimeField(OldTimeField<FieldType>&& otf);


    //- Destructor
    ~OldTimeField();


    // Member Functions

        //- Return the time index of the field
        inline label timeIndex() const;

        //- Return a non-const reference to the time index of the field
        inline label& timeIndex();

        //- Return whether or not this is an old-time field
        bool isOldTime() const;

        //- Store the old-time fields
        void storeOldTimes() const;

        //- Clear old time fields
        void clearOldTimes();

        //- Set the oldest field pointer to nullObjectPtr. This removes the
        //  field whilst maintaining a tag distinct from nullptr so that the
        //  field can be reinstated on the next storeOldTimes.
        void nullOldestTime();

        //- Return the number of old time fields stored
        label nOldTimes(const bool includeNull=true) const;

        //- Return the old time field
        const FieldType& oldTime() const;

        //- Return a non-const reference to the old time field
        FieldType& oldTimeRef();

        //- Return the n-th old time field
        const FieldType& oldTime(const label n) const;

        //- Return a non-const reference to the n-th old time field
        FieldType& oldTimeRef(const label n);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const OldTimeField<FieldType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define USING_OLD_TIME_FIELD(FieldType)                                        \
    using OldTimeField<FieldType>::readOldTimeIfPresent;                       \
    using OldTimeField<FieldType>::copyOldTimes;                               \
    using OldTimeField<FieldType>::timeIndex;                                  \
    using OldTimeField<FieldType>::isOldTime;                                  \
    using OldTimeField<FieldType>::storeOldTimes;                              \
    using OldTimeField<FieldType>::clearOldTimes;                              \
    using OldTimeField<FieldType>::nullOldestTime;                             \
    using OldTimeField<FieldType>::nOldTimes;                                  \
    using OldTimeField<FieldType>::oldTime;                                    \
    using OldTimeField<FieldType>::oldTimeRef;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "OldTimeFieldI.H"

#ifdef NoRepository
    #include "OldTimeField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
