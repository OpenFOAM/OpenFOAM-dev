/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "Pstream.H"
#include "PstreamBuffers.H"
#include "PstreamCombineReduceOps.H"
#include "flipOp.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, class CombineOp, class negateOp>
void Foam::mapDistributeBase::flipAndCombine
(
    const UList<label>& map,
    const bool hasFlip,
    const UList<T>& rhs,
    const CombineOp& cop,
    const negateOp& negOp,
    List<T>& lhs
)
{
    if (hasFlip)
    {
        forAll(map, i)
        {
            if (map[i] > 0)
            {
                label index = map[i]-1;
                cop(lhs[index], rhs[i]);
            }
            else if (map[i] < 0)
            {
                label index = -map[i]-1;
                cop(lhs[index], negOp(rhs[i]));
            }
            else
            {
                FatalErrorInFunction
                    << "At index " << i << " out of " << map.size()
                    << " have illegal index " << map[i]
                    << " for field " << rhs.size() << " with flipMap"
                    << exit(FatalError);
            }
        }
    }
    else
    {
        forAll(map, i)
        {
            cop(lhs[map[i]], rhs[i]);
        }
    }
}


template<class T, class negateOp>
T Foam::mapDistributeBase::accessAndFlip
(
    const UList<T>& fld,
    const label index,
    const bool hasFlip,
    const negateOp& negOp
)
{
    T t;
    if (hasFlip)
    {
        if (index > 0)
        {
            t = fld[index-1];
        }
        else if (index < 0)
        {
            t = negOp(fld[-index-1]);
        }
        else
        {
            FatalErrorInFunction
                << "Illegal index " << index
                << " into field of size " << fld.size()
                << " with face-flipping"
                << exit(FatalError);
            t = fld[index];
        }
    }
    else
    {
        t = fld[index];
    }
    return t;
}


// Distribute list.
template<class T, class negateOp>
void Foam::mapDistributeBase::distribute
(
    const Pstream::commsTypes commsType,
    const List<labelPair>& schedule,
    const label constructSize,
    const labelListList& subMap,
    const bool subHasFlip,
    const labelListList& constructMap,
    const bool constructHasFlip,
    List<T>& field,
    const negateOp& negOp,
    const int tag
)
{
    if (!Pstream::parRun())
    {
        // Do only me to me.

        const labelList& mySubMap = subMap[Pstream::myProcNo()];

        List<T> subField(mySubMap.size());
        forAll(mySubMap, i)
        {
            subField[i] = accessAndFlip(field, mySubMap[i], subHasFlip, negOp);
        }

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        field.setSize(constructSize);

        flipAndCombine
        (
            map,
            constructHasFlip,
            subField,
            eqOp<T>(),
            negOp,
            field
        );

        return;
    }

    if (commsType == Pstream::commsTypes::blocking)
    {
        // Since buffered sending can reuse the field to collect the
        // received data.

        // Send sub field to neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = subMap[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                OPstream toNbr(Pstream::commsTypes::blocking, domain, 0, tag);

                List<T> subField(map.size());
                forAll(subField, i)
                {
                    subField[i] = accessAndFlip
                    (
                        field,
                        map[i],
                        subHasFlip,
                        negOp
                    );
                }
                toNbr << subField;
            }
        }

        // Subset myself
        const labelList& mySubMap = subMap[Pstream::myProcNo()];

        List<T> subField(mySubMap.size());
        forAll(mySubMap, i)
        {
            subField[i] = accessAndFlip(field, mySubMap[i], subHasFlip, negOp);
        }

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        field.setSize(constructSize);

        flipAndCombine
        (
            map,
            constructHasFlip,
            subField,
            eqOp<T>(),
            negOp,
            field
        );

        // Receive sub field from neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = constructMap[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                IPstream fromNbr(Pstream::commsTypes::blocking, domain, 0, tag);
                List<T> subField(fromNbr);

                checkReceivedSize(domain, map.size(), subField.size());

                flipAndCombine
                (
                    map,
                    constructHasFlip,
                    subField,
                    eqOp<T>(),
                    negOp,
                    field
                );
            }
        }
    }
    else if (commsType == Pstream::commsTypes::scheduled)
    {
        // Need to make sure I don't overwrite field with received data
        // since the data might need to be sent to another processor. So
        // allocate a new field for the results.
        List<T> newField(constructSize);

        // Receive sub field from myself
        {
            const labelList& mySubMap = subMap[Pstream::myProcNo()];

            List<T> subField(mySubMap.size());
            forAll(subField, i)
            {
                subField[i] = accessAndFlip
                (
                    field,
                    mySubMap[i],
                    subHasFlip,
                    negOp
                );
            }

            // Receive sub field from myself (subField)
            flipAndCombine
            (
                constructMap[Pstream::myProcNo()],
                constructHasFlip,
                subField,
                eqOp<T>(),
                negOp,
                newField
            );
        }

        // Schedule will already have pruned 0-sized comms
        forAll(schedule, i)
        {
            const labelPair& twoProcs = schedule[i];
            // twoProcs is a swap pair of processors. The first one is the
            // one that needs to send first and then receive.

            label sendProc = twoProcs[0];
            label recvProc = twoProcs[1];

            if (Pstream::myProcNo() == sendProc)
            {
                // I am send first, receive next
                {
                    OPstream toNbr
                    (
                        Pstream::commsTypes::scheduled,
                        recvProc,
                        0,
                        tag
                    );

                    const labelList& map = subMap[recvProc];
                    List<T> subField(map.size());
                    forAll(subField, i)
                    {
                        subField[i] = accessAndFlip
                        (
                            field,
                            map[i],
                            subHasFlip,
                            negOp
                        );
                    }
                    toNbr << subField;
                }
                {
                    IPstream fromNbr
                    (
                        Pstream::commsTypes::scheduled,
                        recvProc,
                        0,
                        tag
                    );
                    List<T> subField(fromNbr);

                    const labelList& map = constructMap[recvProc];

                    checkReceivedSize(recvProc, map.size(), subField.size());

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        subField,
                        eqOp<T>(),
                        negOp,
                        newField
                    );
                }
            }
            else
            {
                // I am receive first, send next
                {
                    IPstream fromNbr
                    (
                        Pstream::commsTypes::scheduled,
                        sendProc,
                        0,
                        tag
                    );
                    List<T> subField(fromNbr);

                    const labelList& map = constructMap[sendProc];

                    checkReceivedSize(sendProc, map.size(), subField.size());

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        subField,
                        eqOp<T>(),
                        negOp,
                        newField
                    );
                }
                {
                    OPstream toNbr
                    (
                        Pstream::commsTypes::scheduled,
                        sendProc,
                        0,
                        tag
                    );

                    const labelList& map = subMap[sendProc];
                    List<T> subField(map.size());
                    forAll(subField, i)
                    {
                        subField[i] = accessAndFlip
                        (
                            field,
                            map[i],
                            subHasFlip,
                            negOp
                        );
                    }
                    toNbr << subField;
                }
            }
        }
        field.transfer(newField);
    }
    else if (commsType == Pstream::commsTypes::nonBlocking)
    {
        label nOutstanding = Pstream::nRequests();

        if (!contiguous<T>())
        {
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking, tag);

            // Stream data into buffer
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    // Put data into send buffer
                    UOPstream toDomain(domain, pBufs);

                    List<T> subField(map.size());
                    forAll(subField, i)
                    {
                        subField[i] = accessAndFlip
                        (
                            field,
                            map[i],
                            subHasFlip,
                            negOp
                        );
                    }
                    toDomain << subField;
                }
            }

            // Start receiving. Do not block.
            pBufs.finishedSends(false);

            {
                // Set up 'send' to myself
                const labelList& mySub = subMap[Pstream::myProcNo()];
                List<T> mySubField(mySub.size());
                forAll(mySub, i)
                {
                    mySubField[i] = accessAndFlip
                    (
                        field,
                        mySub[i],
                        subHasFlip,
                        negOp
                    );
                }
                // Combine bits. Note that can reuse field storage
                field.setSize(constructSize);
                // Receive sub field from myself
                {
                    const labelList& map = constructMap[Pstream::myProcNo()];

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        mySubField,
                        eqOp<T>(),
                        negOp,
                        field
                    );
                }
            }

            // Block ourselves, waiting only for the current comms
            Pstream::waitRequests(nOutstanding);

            // Consume
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    UIPstream str(domain, pBufs);
                    List<T> recvField(str);

                    checkReceivedSize(domain, map.size(), recvField.size());

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        recvField,
                        eqOp<T>(),
                        negOp,
                        field
                    );
                }
            }
        }
        else
        {
            // Set up sends to neighbours

            List<List<T>> sendFields(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    List<T>& subField = sendFields[domain];
                    subField.setSize(map.size());
                    forAll(map, i)
                    {
                        subField[i] = accessAndFlip
                        (
                            field,
                            map[i],
                            subHasFlip,
                            negOp
                        );
                    }

                    OPstream::write
                    (
                        Pstream::commsTypes::nonBlocking,
                        domain,
                        reinterpret_cast<const char*>(subField.begin()),
                        subField.byteSize(),
                        tag
                    );
                }
            }

            // Set up receives from neighbours

            List<List<T>> recvFields(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    recvFields[domain].setSize(map.size());
                    IPstream::read
                    (
                        Pstream::commsTypes::nonBlocking,
                        domain,
                        reinterpret_cast<char*>(recvFields[domain].begin()),
                        recvFields[domain].byteSize(),
                        tag
                    );
                }
            }


            // Set up 'send' to myself

            {
                const labelList& map = subMap[Pstream::myProcNo()];

                List<T>& subField = sendFields[Pstream::myProcNo()];
                subField.setSize(map.size());
                forAll(map, i)
                {
                    subField[i] = accessAndFlip
                    (
                        field,
                        map[i],
                        subHasFlip,
                        negOp
                    );
                }
            }


            // Combine bits. Note that can reuse field storage

            field.setSize(constructSize);


            // Receive sub field from myself (sendFields[Pstream::myProcNo()])
            {
                const labelList& map = constructMap[Pstream::myProcNo()];
                const List<T>& subField = sendFields[Pstream::myProcNo()];

                flipAndCombine
                (
                    map,
                    constructHasFlip,
                    subField,
                    eqOp<T>(),
                    negOp,
                    field
                );
            }


            // Wait for all to finish

            Pstream::waitRequests(nOutstanding);


            // Collect neighbour fields

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    const List<T>& subField = recvFields[domain];

                    checkReceivedSize(domain, map.size(), subField.size());

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        subField,
                        eqOp<T>(),
                        negOp,
                        field
                    );
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown communication schedule " << int(commsType)
            << abort(FatalError);
    }
}


// Distribute list.
template<class T, class CombineOp, class negateOp>
void Foam::mapDistributeBase::distribute
(
    const Pstream::commsTypes commsType,
    const List<labelPair>& schedule,
    const label constructSize,
    const labelListList& subMap,
    const bool subHasFlip,
    const labelListList& constructMap,
    const bool constructHasFlip,
    List<T>& field,
    const CombineOp& cop,
    const negateOp& negOp,
    const T& nullValue,
    const int tag
)
{
    if (!Pstream::parRun())
    {
        // Do only me to me.

        const labelList& mySubMap = subMap[Pstream::myProcNo()];

        List<T> subField(mySubMap.size());
        forAll(mySubMap, i)
        {
            subField[i] = accessAndFlip(field, mySubMap[i], subHasFlip, negOp);
        }

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        field.setSize(constructSize);
        field = nullValue;

        flipAndCombine(map, constructHasFlip, subField, cop, negOp, field);

        return;
    }

    if (commsType == Pstream::commsTypes::blocking)
    {
        // Since buffered sending can reuse the field to collect the
        // received data.

        // Send sub field to neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = subMap[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                OPstream toNbr(Pstream::commsTypes::blocking, domain, 0, tag);
                List<T> subField(map.size());
                forAll(subField, i)
                {
                    subField[i] = accessAndFlip
                    (
                        field,
                        map[i],
                        subHasFlip,
                        negOp
                    );
                }
                toNbr << subField;
            }
        }

        // Subset myself
        const labelList& mySubMap = subMap[Pstream::myProcNo()];

        List<T> subField(mySubMap.size());
        forAll(mySubMap, i)
        {
            subField[i] = accessAndFlip(field, mySubMap[i], subHasFlip, negOp);
        }

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        field.setSize(constructSize);
        field = nullValue;

        flipAndCombine(map, constructHasFlip, subField, cop, negOp, field);

        // Receive sub field from neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = constructMap[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                IPstream fromNbr(Pstream::commsTypes::blocking, domain, 0, tag);
                List<T> subField(fromNbr);

                checkReceivedSize(domain, map.size(), subField.size());

                flipAndCombine
                (
                    map,
                    constructHasFlip,
                    subField,
                    cop,
                    negOp,
                    field
                );
            }
        }
    }
    else if (commsType == Pstream::commsTypes::scheduled)
    {
        // Need to make sure I don't overwrite field with received data
        // since the data might need to be sent to another processor. So
        // allocate a new field for the results.
        List<T> newField(constructSize, nullValue);

        {
            const labelList& mySubMap = subMap[Pstream::myProcNo()];

            // Subset myself
            List<T> subField(mySubMap.size());
            forAll(subField, i)
            {
                subField[i] = accessAndFlip
                (
                    field,
                    mySubMap[i],
                    subHasFlip,
                    negOp
                );
            }

            // Receive sub field from myself (subField)
            const labelList& map = constructMap[Pstream::myProcNo()];

            flipAndCombine
            (
                map,
                constructHasFlip,
                subField,
                cop,
                negOp,
                newField
            );
        }


        // Schedule will already have pruned 0-sized comms
        forAll(schedule, i)
        {
            const labelPair& twoProcs = schedule[i];
            // twoProcs is a swap pair of processors. The first one is the
            // one that needs to send first and then receive.

            label sendProc = twoProcs[0];
            label recvProc = twoProcs[1];

            if (Pstream::myProcNo() == sendProc)
            {
                // I am send first, receive next
                {
                    OPstream toNbr
                    (
                        Pstream::commsTypes::scheduled,
                        recvProc,
                        0,
                        tag
                    );

                    const labelList& map = subMap[recvProc];

                    List<T> subField(map.size());
                    forAll(subField, i)
                    {
                        subField[i] = accessAndFlip
                        (
                            field,
                            map[i],
                            subHasFlip,
                            negOp
                        );
                    }
                    toNbr << subField;
                }
                {
                    IPstream fromNbr
                    (
                        Pstream::commsTypes::scheduled,
                        recvProc,
                        0,
                        tag
                    );
                    List<T> subField(fromNbr);
                    const labelList& map = constructMap[recvProc];

                    checkReceivedSize(recvProc, map.size(), subField.size());

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        subField,
                        cop,
                        negOp,
                        newField
                    );
                }
            }
            else
            {
                // I am receive first, send next
                {
                    IPstream fromNbr
                    (
                        Pstream::commsTypes::scheduled,
                        sendProc,
                        0,
                        tag
                    );
                    List<T> subField(fromNbr);
                    const labelList& map = constructMap[sendProc];

                    checkReceivedSize(sendProc, map.size(), subField.size());

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        subField,
                        cop,
                        negOp,
                        newField
                    );
                }
                {
                    OPstream toNbr
                    (
                        Pstream::commsTypes::scheduled,
                        sendProc,
                        0,
                        tag
                    );

                    const labelList& map = subMap[sendProc];

                    List<T> subField(map.size());
                    forAll(subField, i)
                    {
                        subField[i] = accessAndFlip
                        (
                            field,
                            map[i],
                            subHasFlip,
                            negOp
                        );
                    }
                    toNbr << subField;
                }
            }
        }
        field.transfer(newField);
    }
    else if (commsType == Pstream::commsTypes::nonBlocking)
    {
        label nOutstanding = Pstream::nRequests();

        if (!contiguous<T>())
        {
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking, tag);

            // Stream data into buffer
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    // Put data into send buffer
                    UOPstream toDomain(domain, pBufs);

                    List<T> subField(map.size());
                    forAll(subField, i)
                    {
                        subField[i] = accessAndFlip
                        (
                            field,
                            map[i],
                            subHasFlip,
                            negOp
                        );
                    }
                    toDomain << subField;
                }
            }

            // Start receiving. Do not block.
            pBufs.finishedSends(false);

            {
                // Set up 'send' to myself
                const labelList& myMap = subMap[Pstream::myProcNo()];

                List<T> mySubField(myMap.size());
                forAll(myMap, i)
                {
                    mySubField[i] = accessAndFlip
                    (
                        field,
                        myMap[i],
                        subHasFlip,
                        negOp
                    );
                }

                // Combine bits. Note that can reuse field storage
                field.setSize(constructSize);
                field = nullValue;
                // Receive sub field from myself
                {
                    const labelList& map = constructMap[Pstream::myProcNo()];

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        mySubField,
                        cop,
                        negOp,
                        field
                    );
                }
            }

            // Block ourselves, waiting only for the current comms
            Pstream::waitRequests(nOutstanding);

            // Consume
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    UIPstream str(domain, pBufs);
                    List<T> recvField(str);

                    checkReceivedSize(domain, map.size(), recvField.size());

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        recvField,
                        cop,
                        negOp,
                        field
                    );
                }
            }
        }
        else
        {
            // Set up sends to neighbours

            List<List<T>> sendFields(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    List<T>& subField = sendFields[domain];
                    subField.setSize(map.size());
                    forAll(map, i)
                    {
                        subField[i] = accessAndFlip
                        (
                            field,
                            map[i],
                            subHasFlip,
                            negOp
                        );
                    }

                    OPstream::write
                    (
                        Pstream::commsTypes::nonBlocking,
                        domain,
                        reinterpret_cast<const char*>(subField.begin()),
                        subField.size()*sizeof(T),
                        tag
                    );
                }
            }

            // Set up receives from neighbours

            List<List<T>> recvFields(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    recvFields[domain].setSize(map.size());
                    UIPstream::read
                    (
                        Pstream::commsTypes::nonBlocking,
                        domain,
                        reinterpret_cast<char*>(recvFields[domain].begin()),
                        recvFields[domain].size()*sizeof(T),
                        tag
                    );
                }
            }

            // Set up 'send' to myself

            {
                const labelList& map = subMap[Pstream::myProcNo()];

                List<T>& subField = sendFields[Pstream::myProcNo()];
                subField.setSize(map.size());
                forAll(map, i)
                {
                    subField[i] = accessAndFlip
                    (
                        field,
                        map[i],
                        subHasFlip,
                        negOp
                    );
                }
            }


            // Combine bits. Note that can reuse field storage

            field.setSize(constructSize);
            field = nullValue;

            // Receive sub field from myself (subField)
            {
                const labelList& map = constructMap[Pstream::myProcNo()];
                const List<T>& subField = sendFields[Pstream::myProcNo()];

                flipAndCombine
                (
                    map,
                    constructHasFlip,
                    subField,
                    cop,
                    negOp,
                    field
                );
            }


            // Wait for all to finish

            Pstream::waitRequests(nOutstanding);


            // Collect neighbour fields

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    const List<T>& subField = recvFields[domain];

                    checkReceivedSize(domain, map.size(), subField.size());

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        subField,
                        cop,
                        negOp,
                        field
                    );
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown communication schedule " << int(commsType)
            << abort(FatalError);
    }
}


template<class T>
void Foam::mapDistributeBase::send(PstreamBuffers& pBufs, const List<T>& field)
const
{
    // Stream data into buffer
    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& map = subMap_[domain];

        if (map.size())
        {
            // Put data into send buffer
            UOPstream toDomain(domain, pBufs);

            List<T> subField(map.size());
            forAll(subField, i)
            {
                subField[i] = accessAndFlip
                (
                    field,
                    map[i],
                    subHasFlip_,
                    flipOp()
                );
            }
            toDomain << subField;
        }
    }

    // Start sending and receiving but do not block.
    pBufs.finishedSends(false);
}


template<class T>
void Foam::mapDistributeBase::receive(PstreamBuffers& pBufs, List<T>& field)
const
{
    // Consume
    field.setSize(constructSize_);

    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& map = constructMap_[domain];

        if (map.size())
        {
            UIPstream str(domain, pBufs);
            List<T> recvField(str);

            if (recvField.size() != map.size())
            {
                FatalErrorInFunction
                    << "Expected from processor " << domain
                    << " " << map.size() << " but received "
                    << recvField.size() << " elements."
                    << abort(FatalError);
            }

            flipAndCombine
            (
                map,
                constructHasFlip_,
                recvField,
                eqOp<T>(),
                flipOp(),
                field
            );
        }
    }
}


//- Distribute data using default commsType.
template<class T, class negateOp>
void Foam::mapDistributeBase::distribute
(
    List<T>& fld,
    const negateOp& negOp,
    const int tag
) const
{
    if (Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking)
    {
        distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            constructSize_,
            subMap_,
            subHasFlip_,
            constructMap_,
            constructHasFlip_,
            fld,
            negOp,
            tag
        );
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        distribute
        (
            Pstream::commsTypes::scheduled,
            schedule(),
            constructSize_,
            subMap_,
            subHasFlip_,
            constructMap_,
            constructHasFlip_,
            fld,
            negOp,
            tag
        );
    }
    else
    {
        distribute
        (
            Pstream::commsTypes::blocking,
            List<labelPair>(),
            constructSize_,
            subMap_,
            subHasFlip_,
            constructMap_,
            constructHasFlip_,
            fld,
            negOp,
            tag
        );
    }
}


//- Distribute data using default commsType.
template<class T>
void Foam::mapDistributeBase::distribute
(
    List<T>& fld,
    const int tag
) const
{
    distribute(fld, flipOp(), tag);
}


//- Distribute data using default commsType.
template<class T>
void Foam::mapDistributeBase::distribute
(
    DynamicList<T>& fld,
    const int tag
) const
{
    fld.shrink();

    List<T>& fldList = static_cast<List<T>& >(fld);

    distribute(fldList, tag);

    fld.setCapacity(fldList.size());
}


//- Reverse distribute data using default commsType.
template<class T>
void Foam::mapDistributeBase::reverseDistribute
(
    const label constructSize,
    List<T>& fld,
    const int tag
) const
{
    if (Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking)
    {
        distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            constructSize,
            constructMap_,
            constructHasFlip_,
            subMap_,
            subHasFlip_,
            fld,
            flipOp(),
            tag
        );
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        distribute
        (
            Pstream::commsTypes::scheduled,
            schedule(),
            constructSize,
            constructMap_,
            constructHasFlip_,
            subMap_,
            subHasFlip_,
            fld,
            flipOp(),
            tag
        );
    }
    else
    {
        distribute
        (
            Pstream::commsTypes::blocking,
            List<labelPair>(),
            constructSize,
            constructMap_,
            constructHasFlip_,
            subMap_,
            subHasFlip_,
            fld,
            flipOp(),
            tag
        );
    }
}


//- Reverse distribute data using default commsType.
//  Since constructSize might be larger than supplied size supply
//  a nullValue
template<class T>
void Foam::mapDistributeBase::reverseDistribute
(
    const label constructSize,
    const T& nullValue,
    List<T>& fld,
    const int tag
) const
{
    if (Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking)
    {
        distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            constructSize,
            constructMap_,
            constructHasFlip_,
            subMap_,
            subHasFlip_,
            fld,
            eqOp<T>(),
            flipOp(),
            nullValue,
            tag
        );
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        distribute
        (
            Pstream::commsTypes::scheduled,
            schedule(),
            constructSize,
            constructMap_,
            constructHasFlip_,
            subMap_,
            subHasFlip_,
            fld,
            eqOp<T>(),
            flipOp(),
            nullValue,
            tag
        );
    }
    else
    {
        distribute
        (
            Pstream::commsTypes::blocking,
            List<labelPair>(),
            constructSize,
            constructMap_,
            constructHasFlip_,
            subMap_,
            subHasFlip_,
            fld,
            eqOp<T>(),
            flipOp(),
            nullValue,
            tag
        );
    }
}


// ************************************************************************* //
