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

#include "Histogram.H"
#include "ListOps.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class List>
void Foam::Histogram<List>::count(const List& bins, const List& l)
{
    if (bins.size() < 2)
    {
        FatalErrorIn("Histogram<List>::count(const List&, const List&)")
            << "Should have at least two values in bins. Now:" << bins
            << exit(FatalError);
    }

    counts_.setSize(bins.size()-1);
    counts_ = 0;

    nLow_ = 0;
    nHigh_ = 0;

    forAll(l, i)
    {
        label index = findLower(bins, l[i]);

        if (index == -1)
        {
            nLow_++;
        }
        else if (index == bins.size()-1)
        {
            nHigh_++;
        }
        else
        {
            counts_[index]++;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class List>
Foam::Histogram<List>::Histogram(const List& bins, const List& l)
{
    count(bins, l);
}


template<class List>
Foam::Histogram<List>::Histogram
(
    const typename List::const_reference min,
    const typename List::const_reference max,
    const label nBins,
    const List& l
)
{
    List bins(nBins+1);

    typename List::value_type span = (max-min) / nBins;

    bins[0] = min;

    for (label i = 1; i < nBins; i++)
    {
        bins[i] = bins[i-1] + span;
    }

    // Set max directly to avoid truncation errors.
    bins[nBins] = max;

    count(bins, l);
}


// ************************************************************************* //
