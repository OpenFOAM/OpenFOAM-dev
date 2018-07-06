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

\*---------------------------------------------------------------------------*/

#include "PrintTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class KeyType, class DataType>
Foam::PrintTable<KeyType, DataType>::PrintTable()
:
    table_(),
    title_(string::null)
{}


template<class KeyType, class DataType>
Foam::PrintTable<KeyType, DataType>::PrintTable(const string& title)
:
    table_(),
    title_(title)
{}


template<class KeyType, class DataType>
Foam::PrintTable<KeyType, DataType>::PrintTable
(
    const PrintTable<KeyType, DataType>& table
)
:
    table_(table.table_),
    title_(table.title_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class KeyType, class DataType>
Foam::PrintTable<KeyType, DataType>::~PrintTable()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class KeyType, class DataType>
void Foam::PrintTable<KeyType, DataType>::print
(
    Ostream& os,
    const bool printSum,
    const bool printAverage
) const
{
    HashTable<HashTable<DataType, label>, KeyType> combinedTable;

    List<HashTableData> procData(Pstream::nProcs(), HashTableData());

    procData[Pstream::myProcNo()] = table_;

    Pstream::gatherList(procData);

    if (Pstream::master())
    {
        label largestKeyLength = 6;
        label largestDataLength = 0;

        List<label> largestProcSize(Pstream::nProcs(), 0);

        forAll(procData, proci)
        {
            const HashTableData& procIData = procData[proci];

            for
            (
                typename HashTableData::const_iterator iter = procIData.begin();
                iter != procIData.end();
                ++iter
            )
            {
                if (!combinedTable.found(iter.key()))
                {
                    combinedTable.insert
                    (
                        iter.key(),
                        HashTable<DataType, label>()
                    );
                }

                HashTable<DataType, label>& key = combinedTable[iter.key()];

                key.insert(proci, iter());

                for
                (
                    typename HashTable<DataType, label>
                        ::const_iterator dataIter = key.begin();
                    dataIter != key.end();
                    ++dataIter
                )
                {
                    std::ostringstream buf;
                    buf << dataIter();

                    largestDataLength = max
                    (
                        largestDataLength,
                        label(buf.str().length())
                    );
                }

                std::ostringstream buf;
                buf << iter.key();

                largestKeyLength = max
                (
                    largestKeyLength,
                    label(buf.str().length())
                );
            }
        }

        os.width(largestKeyLength);
        os  << nl << indent << tab << "# " << title_.c_str() << endl;

        os.width(largestKeyLength);
        os  << indent << "# Proc";

        forAll(procData, proci)
        {
            os  << tab;
            os.width(largestDataLength);
            os  << proci;
        }

        if (printSum)
        {
            os  << tab;
            os.width(largestDataLength);
            os  << "Sum";
        }

        if (printAverage)
        {
            os  << tab;
            os.width(largestDataLength);
            os  << "Average";
        }

        os  << endl;

        const List<KeyType>& sortedTable = combinedTable.sortedToc();

        forAll(sortedTable, keyI)
        {
            const HashTable<DataType, label>& procDataList
                = combinedTable[sortedTable[keyI]];

            os.width(largestKeyLength);
            os  << indent << sortedTable[keyI];

            forAll(procDataList, elemI)
            {
                os  << tab;
                os.width(largestDataLength);
                os  << procDataList[elemI];
            }

            if (printSum)
            {
                DataType sum = 0;
                forAll(procDataList, elemI)
                {
                    sum += procDataList[elemI];
                }

                os  << tab;
                os.width(largestDataLength);
                os  << sum;

                if (printAverage)
                {
                    os  << tab;
                    os.width(largestDataLength);
                    os  << sum/Pstream::nProcs();
                }
            }

            os  << endl;
        }
    }
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class KeyType, class DataType>
void Foam::PrintTable<KeyType, DataType>::operator=
(
    const PrintTable<KeyType, DataType>& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


// ************************************************************************* //
