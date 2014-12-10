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

Class
    Foam::Distribution

Description
    Accumulating histogram of component values.
    Specified bin resolution, automatic generation of bins.

SourceFiles
    DistributionI.H
    Distribution.C
    DistributionIO.C

\*---------------------------------------------------------------------------*/

#ifndef Distribution_H
#define Distribution_H

#include "List.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class Distribution;

template<class Type>
Istream& operator>>(Istream&, Distribution<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const Distribution<Type>&);

/*---------------------------------------------------------------------------*\
                         Class Distribution Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Distribution
:
    public List< List<scalar> >
{
    // Private data

        //- Width of the bin for each component
        Type binWidth_;

        //- The start bin index of each component
        List<label> listStarts_;


public:

    //- Component type
    typedef typename pTraits<Type>::cmptType cmptType;


    // Constructors

        //- Construct null
        Distribution();

        //- Construct from separate binWidth for each component
        Distribution(const Type& binWidth);

        //- Construct as copy
        Distribution(const Distribution& d);


    //- Destructor
        ~Distribution();

    // Member Functions

        //- Sum the total weight added to the component in the
        //  argument
        scalar totalWeight(direction cmpt) const;

        List<label> keys(direction cmpt) const;

        //- Return the appropriate List index for the given bin index.
        //  Resizes the List if required
        label index(direction cmpt, label n);

        //- Returns the indices of the first and last non-zero entries
        Pair<label> validLimits(direction cmpt) const;

        Type mean() const;

        // From http://mathworld.wolfram.com/StatisticalMedian.html
        // The statistical median is the value of the Distribution
        // variable where the cumulative Distribution = 0.5.
        Type median() const;

        //- Add a value to the distribution, optionally specifying a weight
        void add
        (
            const Type& valueToAdd,
            const Type& weight = pTraits<Type>::one
        );

        //- Return the normalised distribution (probability density)
        //  and bins
        List< List<Pair<scalar> > > normalised() const;

        //- Return the distribution of the total bin weights
        List< List < Pair<scalar> > > raw() const;

        //- Return the cumulative normalised distribution and
        //  integration locations (at end of bins)
        List< List<Pair<scalar> > > cumulativeNormalised() const;

        //- Return the cumulative total bin weights and integration
        //  locations (at end of bins)
        List< List<Pair<scalar> > > cumulativeRaw() const;

        //- Resets the Distribution by clearing the stored lists.
        //  Leaves the same number of them and the same binWidth.
        void clear();


    // Access

        //- Return the bin width
        inline const Type& binWidth() const;

        //- Return the List start bin indices
        inline const List<label>& listStarts() const;

    // Write

        //- Write the distribution to file: key normalised raw.
        //  Produces a separate file for each component.
        void write(const fileName& filePrefix) const;


    // Member Operators

        void operator=(const Distribution<Type>&);

    // IOstream Operators

        friend Istream& operator>> <Type>
        (
            Istream&,
            Distribution<Type>&
        );

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const Distribution<Type>&
        );
};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Type>
Distribution<Type> operator+
(
    const Distribution<Type>&,
    const Distribution<Type>&
);


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DistributionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Distribution.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
