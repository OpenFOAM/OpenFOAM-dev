/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::Histogram

Description
    Calculates the counts per bin of a list.

SourceFiles
    Histogram.C

\*---------------------------------------------------------------------------*/

#ifndef Histogram_H
#define Histogram_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                          Class Histogram Declaration
\*---------------------------------------------------------------------------*/

template<class List>
class Histogram
{
    // Private data

        //- Counts per bin
        labelList counts_;

        //- Number of <= lowest bin
        label nLow_;

        //- Number of > highest bin
        label nHigh_;


    // Private Member Functions

        void count(const List& bins, const List& l);

        //- Disallow default bitwise copy construct
        Histogram(const Histogram&);

        //- Disallow default bitwise assignment
        void operator=(const Histogram&);


public:

    // Constructors

        //- Construct given bin values and input list
        Histogram(const List& bins, const List& l);

        //- Construct given min, max, number of bins and input list
        Histogram
        (
            const typename List::const_reference min,
            const typename List::const_reference max,
            const label nBins,
            const List& l
        );


        // Access

            //- Return the counts per bin
            inline const labelList& counts() const
            {
                return counts_;
            }

            //- Return the number of elements <= bins[0]
            //  (so inclusive lowest bin value)
            inline label nLow() const
            {
                return nLow_;
            }

            //- Return the number of elements > bins[bins.size()-1]
            //  (so exclusive highest bin value)
            inline label nHigh() const
            {
                return nHigh_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//    #include "HistogramI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Histogram.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
