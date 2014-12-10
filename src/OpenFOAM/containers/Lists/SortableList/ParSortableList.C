/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "ParSortableList.H"
#include "SortableList.H"
#include "Pstream.H"
#include "ListListOps.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::ParSortableList<Type>::write
(
    const List<Type>& elems,
    Ostream& os
) const
{
    os << '(';

    forAll(elems, elemI)
    {
        os << ' ' << elems[elemI];
    }
    os << ')';
}


// Copy src, starting at destI into dest.
template<class Type>
void Foam::ParSortableList<Type>::copyInto
(
    const List<Type>& values,
    const labelList& indices,
    const label fromProcNo,
    label& destI,
    List<taggedValue>& dest
) const
{
    forAll(values, elemI)
    {
        taggedValue& tagVal = dest[destI];

        tagVal.value() = values[elemI];
        tagVal.index() = indices[elemI];
        tagVal.procID() = fromProcNo;

        destI++;
    }
}


template<class Type>
void Foam::ParSortableList<Type>::getPivots
(
    const List<Type>& elems,
    List<Type>& pivots
) const
{
    pivots.setSize(Pstream::nProcs());

    label pivotPos = 0;

    forAll(pivots, pivotI)
    {
        pivots[pivotI] = elems[pivotPos];

        pivotPos += elems.size()/Pstream::nProcs();
    }
}


template<class Type>
void Foam::ParSortableList<Type>::checkAndSend
(
    List<Type>& values,
    labelList& indices,
    const label bufSize,
    const label destProcI
) const
{
    if (destProcI != Pstream::myProcNo())
    {
        values.setSize(bufSize);
        indices.setSize(bufSize);

        if (debug)
        {
            Pout<<  "Sending to " << destProcI << " elements:" << values
                << endl;
        }

        {
            OPstream toSlave(Pstream::blocking, destProcI);
            toSlave << values << indices;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from List, sorting the elements
template<class Type>
Foam::ParSortableList<Type>::ParSortableList(const UList<Type>& values)
:
    List<Type>(values),
    indices_(0),
    procs_(0)
{
    sort();
}


// Construct given size. Sort later on.
template<class Type>
Foam::ParSortableList<Type>::ParSortableList(const label size)
:
    List<Type>(size),
    indices_(0),
    procs_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Sort
template<class Type>
void Foam::ParSortableList<Type>::sort()
{
    //
    // 0. Get total size of dataset.
    //

    label n = this->size();

    reduce(n, sumOp<label>());


    // 1. Sort list locally
    SortableList<Type> sorted(*this);

    // Collect elements at pivot points
    labelListList sortedGatherList(Pstream::nProcs());

    labelList& pivots = sortedGatherList[Pstream::myProcNo()];

    getPivots(sorted, pivots);

    if (debug)
    {
        Pout<< "pivots:";
        write(pivots, Pout);
        Pout<< endl;
    }


    //
    // 2. Combine pivotlist per processor onto master, sort, get pivots.
    //

    Pstream::gatherList(sortedGatherList);

    if (Pstream::master())
    {
        labelList allPivots =
            ListListOps::combine<labelList>
            (
                sortedGatherList,
                accessOp<labelList>()
            );

        SortableList<Type> sortedPivots(allPivots);

        if (debug)
        {
            Pout<< "allPivots:";
            write(allPivots, Pout);
            Pout<< endl;
        }

        getPivots(sortedPivots, pivots);
    }
    Pstream::scatter(pivots);

    if (debug)
    {
        Pout<< "new pivots:";
        write(pivots, Pout);
        Pout<< endl;
    }


    //
    // 3. Distribute pivots & distribute.
    //

    label pivotI = 1;
    label destProcI = 0;

    // Buffer for my own data. Keep original index together with value.
    labelList ownValues(sorted.size());
    labelList ownIndices(sorted.size());
    label ownI = 0;

    // Buffer for sending data
    labelList sendValues(sorted.size());
    labelList sendIndices(sorted.size());
    label sendI = 0;

    forAll(sorted, sortedI)
    {
        if ((pivotI < Pstream::nProcs()) && (sorted[sortedI] > pivots[pivotI]))
        {
            checkAndSend(sendValues, sendIndices, sendI, destProcI);

            // Reset buffer.
            sendValues.setSize(sorted.size());
            sendIndices.setSize(sorted.size());
            sendI = 0;

            pivotI++;
            destProcI++;
        }

        if (destProcI != Pstream::myProcNo())
        {
            sendValues[sendI] = sorted[sortedI];
            sendIndices[sendI] = sorted.indices()[sortedI];
            sendI++;
        }
        else
        {
            ownValues[ownI] = sorted[sortedI];
            ownIndices[ownI] = sorted.indices()[sortedI];
            ownI++;
        }
    }


    // Handle trailing send buffer
    if (sendI != 0)
    {
        checkAndSend(sendValues, sendIndices, sendI, destProcI);
    }

    // Print ownValues
    ownValues.setSize(ownI);
    ownIndices.setSize(ownI);

    if (debug & 2)
    {
        Pout<< "Not sending (to myself) elements "
            << ownValues << endl;
    }

    //
    // 4. Combine pieces from all processors & sort. Use indices() from
    // SortableList to remember source processor number.
    //

    // Allocate receive buffer. Acc. to paper upper bound is 2*n/p
    // (n=total size, p=nProcs). Resize later on.
    List<taggedValue> combinedValues(2 * n/Pstream::nProcs());

    label combinedI = 0;

    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        if (procI == Pstream::myProcNo())
        {
            if (debug & 2)
            {
                Pout<< "Copying from own:" << ownValues << endl;
            }

            // Copy ownValues,ownIndices into combined buffer
            copyInto(ownValues, ownIndices, procI, combinedI, combinedValues);
        }
        else
        {
            labelList recValues;
            labelList recIndices;

            {
                if (debug)
                {
                    Pout<< "Receiving from " << procI << endl;
                }

                IPstream fromSlave(Pstream::blocking, procI);

                fromSlave >> recValues >> recIndices;

                if (debug & 2)
                {
                    Pout<< "Received from " << procI
                        << " elements:" << recValues << endl;
                }
            }

            if (debug)
            {
                Pout<< "Copying starting at:" << combinedI << endl;
            }
            copyInto(recValues, recIndices, procI, combinedI, combinedValues);
        }
    }
    combinedValues.setSize(combinedI);

    // Sort according to values
    Foam::sort(combinedValues);

    // Copy into *this
    this->setSize(combinedI);
    indices_.setSize(combinedI);
    procs_.setSize(combinedI);

    forAll(combinedValues, elemI)
    {
        this->operator[](elemI) = combinedValues[elemI].value();
        indices_[elemI] = combinedValues[elemI].index();
        procs_[elemI] = combinedValues[elemI].procID();
    }
}


// ************************************************************************* //
