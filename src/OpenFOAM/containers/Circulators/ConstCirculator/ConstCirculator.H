/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    Foam::ConstCirculator

Description
    Walks over a container as if it were circular. The container must have the
    following members defined:
        - value_type
        - size_type
        - difference_type
        - const_iterator
        - const_reference

    Examples:

    \code
        face f(identity(5));

        // Construct circulator from the face
        ConstCirculator<face> circ(f);

        // First check that the circulator has a size to iterate over.
        // Then circulate around the list starting and finishing at the fulcrum.
        if (circ.size()) do
        {
            Info<< "Iterate forwards over face : " << circ() << endl;

        } while (circ.circulate(CirculatorBase::CLOCKWISE));
    \endcode

    \code
        face f(identity(5));

        ConstCirculator<face> circClockwise(f);
        ConstCirculator<face> circAnticlockwise(f);

        if (circClockwise.size() && circAnticlockwise.size()) do
        {
            Info<< "Iterate forward over face :" << circClockwise() << endl;
            Info<< "Iterate backward over face:" << circAnticlockwise() << endl;
        }
        while
        (
            circClockwise.circulate(CirculatorBase::CLOCKWISE),
            circAnticlockwise.circulate(CirculatorBase::ANTICLOCKWISE)
        );
    \endcode

SourceFiles
    ConstCirculatorI.H

\*---------------------------------------------------------------------------*/

#ifndef ConstCirculator_H
#define ConstCirculator_H

#include "CirculatorBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                      Class ConstCirculator Declaration
\*---------------------------------------------------------------------------*/

template<class ContainerType>
class ConstCirculator
:
    public CirculatorBase
{

protected:

    // Protected data

        //- Iterator pointing to the beginning of the container
        typename ContainerType::const_iterator begin_;

        //- Iterator pointing to the end of the container
        typename ContainerType::const_iterator end_;

        //- Iterator
        typename ContainerType::const_iterator iter_;

        //- Iterator holding the location of the fulcrum (start and end) of
        //  the container. Used to decide when the iterator should stop
        //  circulating over the container
        typename ContainerType::const_iterator fulcrum_;


public:

    // STL type definitions

        //- Type of values ContainerType contains.
        typedef typename ContainerType::value_type      value_type;

        //- The type that can represent the size of ContainerType
        typedef typename ContainerType::size_type       size_type;

        //- The type that can represent the difference between any two
        //  iterator objects.
        typedef typename ContainerType::difference_type difference_type;

        //- Random access iterator for traversing ContainerType.
        typedef typename ContainerType::const_iterator  const_iterator;

        //- Type that can be used for storing into
        //  const ContainerType::value_type objects.
        typedef typename ContainerType::const_reference const_reference;


    // Constructors

        //- Construct null
        inline ConstCirculator();

        //- Construct from a container.
        inline explicit ConstCirculator(const ContainerType& container);

        //- Construct from two iterators
        inline ConstCirculator
        (
            const const_iterator& begin,
            const const_iterator& end
        );

        //- Construct as copy
        inline ConstCirculator(const ConstCirculator<ContainerType>&);


    //- Destructor
    ~ConstCirculator();


    // Member Functions

        //- Return the range of the iterator
        inline size_type size() const;

        //- Circulate around the list in the given direction
        inline bool circulate(const CirculatorBase::direction dir = NONE);

        //- Set the fulcrum to the current position of the iterator
        inline void setFulcrumToIterator();

        //- Set the iterator to the current position of the fulcrum
        inline void setIteratorToFulcrum();

        //- Return the distance between the iterator and the fulcrum. This is
        //  equivalent to the number of rotations of the circulator.
        inline difference_type nRotations() const;

        //- Dereference the next iterator and return
        inline const_reference next() const;

        //- Dereference the previous iterator and return
        inline const_reference prev() const;


    // Member Operators

        //- Assignment operator for circulators that operate on the same
        //  container type
        inline void operator=(const ConstCirculator<ContainerType>&);

        //- Prefix increment. Increments the iterator.
        //  Sets the iterator to the beginning of the container if it reaches
        //  the end
        inline ConstCirculator<ContainerType>& operator++();

        //- Postfix increment. Increments the iterator.
        //  Sets the iterator to the beginning of the container if it reaches
        //  the end
        inline ConstCirculator<ContainerType> operator++(int);

        //- Prefix decrement. Decrements the iterator.
        //  Sets the iterator to the end of the container if it reaches
        //  the beginning
        inline ConstCirculator<ContainerType>& operator--();

        //- Postfix decrement. Decrements the iterator.
        //  Sets the iterator to the end of the container if it reaches
        //  the beginning
        inline ConstCirculator<ContainerType> operator--(int);

        //- Check for equality of this iterator with another iterator that
        //  operate on the same container type
        inline bool operator==(const ConstCirculator<ContainerType>& c) const;

        //- Check for inequality of this iterator with another iterator that
        //  operate on the same container type
        inline bool operator!=(const ConstCirculator<ContainerType>& c) const;

        //- Dereference the iterator and return
        inline const_reference operator*() const;

        //- Dereference the iterator and return
        inline const_reference operator()() const;

        //- Return the difference between this iterator and another iterator
        //  that operate on the same container type
        inline difference_type operator-
        (
            const ConstCirculator<ContainerType>& c
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ConstCirculatorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
