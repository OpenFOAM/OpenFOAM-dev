/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Foam::PrintTable

Description
    Print a table in parallel, e.g.;

    \verbatim
                            Vertex Type Information
                   Proc #        0       1       2       3
                    Total   145680  145278  145751  145359
               Unassigned        0       0       0       0
     nExternalFeatureEdge      883     829     828     960
    nExternalFeaturePoint        8      10      10      12
         nExternalSurface     9533    9488    9502    9482
                     nFar        0       0       0       0
                nInternal   125494  125198  125642  125174
     nInternalFeatureEdge      238     241     241     240
    nInternalFeaturePoint        2       2       2       2
    nInternalNearBoundary        0       0       0       0
         nInternalSurface     9522    9510    9526    9489
                nReferred     7545    7497    7500    7587
    \endverbatim

SourceFiles
    PrintTableI.H
    PrintTable.C

\*---------------------------------------------------------------------------*/

#ifndef PrintTable_H
#define PrintTable_H

#include "HashTable.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class PrintTable Declaration
\*---------------------------------------------------------------------------*/

template<class KeyType, class DataType>
class PrintTable
{
    typedef HashTable<DataType, KeyType> HashTableData;


    // Private data

        //- Hash table holding the data
        HashTableData table_;

        //- Title of the table
        string title_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const PrintTable<KeyType, DataType>&);


public:

    // Constructors

        //- Null constructor
        PrintTable();

        //- Construct with a title
        explicit PrintTable(const string& title);

        //- Copy constructor
        PrintTable(const PrintTable<KeyType, DataType>& table);


    //- Destructor
    ~PrintTable();


    // Member Functions

        //- Add an entry (D) to the given key(K)
        void add(const KeyType& K, const DataType& D);

        //- Print the table
        void print
        (
            Ostream& os,
            const bool printSum = false,
            const bool printAverage = false
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PrintTableI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PrintTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
