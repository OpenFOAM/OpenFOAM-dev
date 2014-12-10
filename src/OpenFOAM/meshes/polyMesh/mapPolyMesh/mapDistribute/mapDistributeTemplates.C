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

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "PstreamBuffers.H"
#include "PstreamCombineReduceOps.H"
#include "globalIndexAndTransform.H"
#include "transformField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Distribute list.
template<class T>
void Foam::mapDistribute::distribute
(
    const Pstream::commsTypes commsType,
    const List<labelPair>& schedule,
    const label constructSize,
    const labelListList& subMap,
    const labelListList& constructMap,
    List<T>& field,
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
            subField[i] = field[mySubMap[i]];
        }

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        field.setSize(constructSize);

        forAll(map, i)
        {
            field[map[i]] = subField[i];
        }
        return;
    }

    if (commsType == Pstream::blocking)
    {
        // Since buffered sending can reuse the field to collect the
        // received data.

        // Send sub field to neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = subMap[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                OPstream toNbr(Pstream::blocking, domain, 0, tag);
                toNbr << UIndirectList<T>(field, map);
            }
        }

        // Subset myself
        const labelList& mySubMap = subMap[Pstream::myProcNo()];

        List<T> subField(mySubMap.size());
        forAll(mySubMap, i)
        {
            subField[i] = field[mySubMap[i]];
        }

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        field.setSize(constructSize);

        forAll(map, i)
        {
            field[map[i]] = subField[i];
        }

        // Receive sub field from neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = constructMap[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                IPstream fromNbr(Pstream::blocking, domain, 0, tag);
                List<T> subField(fromNbr);

                checkReceivedSize(domain, map.size(), subField.size());

                forAll(map, i)
                {
                    field[map[i]] = subField[i];
                }
            }
        }
    }
    else if (commsType == Pstream::scheduled)
    {
        // Need to make sure I don't overwrite field with received data
        // since the data might need to be sent to another processor. So
        // allocate a new field for the results.
        List<T> newField(constructSize);

        // Subset myself
        UIndirectList<T> subField(field, subMap[Pstream::myProcNo()]);

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        forAll(map, i)
        {
            newField[map[i]] = subField[i];
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
                    OPstream toNbr(Pstream::scheduled, recvProc, 0, tag);
                    toNbr << UIndirectList<T>(field, subMap[recvProc]);
                }
                {
                    IPstream fromNbr(Pstream::scheduled, recvProc, 0, tag);
                    List<T> subField(fromNbr);

                    const labelList& map = constructMap[recvProc];

                    checkReceivedSize(recvProc, map.size(), subField.size());

                    forAll(map, i)
                    {
                        newField[map[i]] = subField[i];
                    }
                }
            }
            else
            {
                // I am receive first, send next
                {
                    IPstream fromNbr(Pstream::scheduled, sendProc, 0, tag);
                    List<T> subField(fromNbr);

                    const labelList& map = constructMap[sendProc];

                    checkReceivedSize(sendProc, map.size(), subField.size());

                    forAll(map, i)
                    {
                        newField[map[i]] = subField[i];
                    }
                }
                {
                    OPstream toNbr(Pstream::scheduled, sendProc, 0, tag);
                    toNbr << UIndirectList<T>(field, subMap[sendProc]);
                }
            }
        }
        field.transfer(newField);
    }
    else if (commsType == Pstream::nonBlocking)
    {
        label nOutstanding = Pstream::nRequests();

        if (!contiguous<T>())
        {
            PstreamBuffers pBufs(Pstream::nonBlocking, tag);

            // Stream data into buffer
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    // Put data into send buffer
                    UOPstream toDomain(domain, pBufs);
                    toDomain << UIndirectList<T>(field, map);
                }
            }

            // Start receiving. Do not block.
            pBufs.finishedSends(false);

            {
                // Set up 'send' to myself
                const labelList& mySubMap = subMap[Pstream::myProcNo()];
                List<T> mySubField(mySubMap.size());
                forAll(mySubMap, i)
                {
                    mySubField[i] = field[mySubMap[i]];
                }
                // Combine bits. Note that can reuse field storage
                field.setSize(constructSize);
                // Receive sub field from myself
                {
                    const labelList& map = constructMap[Pstream::myProcNo()];

                    forAll(map, i)
                    {
                        field[map[i]] = mySubField[i];
                    }
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

                    forAll(map, i)
                    {
                        field[map[i]] = recvField[i];
                    }
                }
            }
        }
        else
        {
            // Set up sends to neighbours

            List<List<T > > sendFields(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    List<T>& subField = sendFields[domain];
                    subField.setSize(map.size());
                    forAll(map, i)
                    {
                        subField[i] = field[map[i]];
                    }

                    OPstream::write
                    (
                        Pstream::nonBlocking,
                        domain,
                        reinterpret_cast<const char*>(subField.begin()),
                        subField.byteSize(),
                        tag
                    );
                }
            }

            // Set up receives from neighbours

            List<List<T > > recvFields(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    recvFields[domain].setSize(map.size());
                    IPstream::read
                    (
                        Pstream::nonBlocking,
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
                    subField[i] = field[map[i]];
                }
            }


            // Combine bits. Note that can reuse field storage

            field.setSize(constructSize);


            // Receive sub field from myself (sendFields[Pstream::myProcNo()])
            {
                const labelList& map = constructMap[Pstream::myProcNo()];
                const List<T>& subField = sendFields[Pstream::myProcNo()];

                forAll(map, i)
                {
                    field[map[i]] = subField[i];
                }
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

                    forAll(map, i)
                    {
                        field[map[i]] = subField[i];
                    }
                }
            }
        }
    }
    else
    {
        FatalErrorIn("mapDistribute::distribute(..)")
            << "Unknown communication schedule " << commsType
            << abort(FatalError);
    }
}


// Distribute list.
template<class T, class CombineOp>
void Foam::mapDistribute::distribute
(
    const Pstream::commsTypes commsType,
    const List<labelPair>& schedule,
    const label constructSize,
    const labelListList& subMap,
    const labelListList& constructMap,
    List<T>& field,
    const CombineOp& cop,
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
            subField[i] = field[mySubMap[i]];
        }

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        field.setSize(constructSize);
        field = nullValue;

        forAll(map, i)
        {
            cop(field[map[i]], subField[i]);
        }
        return;
    }

    if (commsType == Pstream::blocking)
    {
        // Since buffered sending can reuse the field to collect the
        // received data.

        // Send sub field to neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = subMap[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                OPstream toNbr(Pstream::blocking, domain, 0, tag);
                toNbr << UIndirectList<T>(field, map);
            }
        }

        // Subset myself
        const labelList& mySubMap = subMap[Pstream::myProcNo()];

        List<T> subField(mySubMap.size());
        forAll(mySubMap, i)
        {
            subField[i] = field[mySubMap[i]];
        }

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        field.setSize(constructSize);
        field = nullValue;

        forAll(map, i)
        {
            cop(field[map[i]], subField[i]);
        }

        // Receive sub field from neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = constructMap[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                IPstream fromNbr(Pstream::blocking, domain, 0, tag);
                List<T> subField(fromNbr);

                checkReceivedSize(domain, map.size(), subField.size());

                forAll(map, i)
                {
                    cop(field[map[i]], subField[i]);
                }
            }
        }
    }
    else if (commsType == Pstream::scheduled)
    {
        // Need to make sure I don't overwrite field with received data
        // since the data might need to be sent to another processor. So
        // allocate a new field for the results.
        List<T> newField(constructSize, nullValue);

        // Subset myself
        UIndirectList<T> subField(field, subMap[Pstream::myProcNo()]);

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        forAll(map, i)
        {
            cop(newField[map[i]], subField[i]);
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
                    OPstream toNbr(Pstream::scheduled, recvProc, 0, tag);
                    toNbr << UIndirectList<T>(field, subMap[recvProc]);
                }
                {
                    IPstream fromNbr(Pstream::scheduled, recvProc, 0, tag);
                    List<T> subField(fromNbr);
                    const labelList& map = constructMap[recvProc];

                    checkReceivedSize(recvProc, map.size(), subField.size());

                    forAll(map, i)
                    {
                        cop(newField[map[i]], subField[i]);
                    }
                }
            }
            else
            {
                // I am receive first, send next
                {
                    IPstream fromNbr(Pstream::scheduled, sendProc, 0, tag);
                    List<T> subField(fromNbr);
                    const labelList& map = constructMap[sendProc];

                    checkReceivedSize(sendProc, map.size(), subField.size());

                    forAll(map, i)
                    {
                        cop(newField[map[i]], subField[i]);
                    }
                }
                {
                    OPstream toNbr(Pstream::scheduled, sendProc, 0, tag);
                    toNbr << UIndirectList<T>(field, subMap[sendProc]);
                }
            }
        }
        field.transfer(newField);
    }
    else if (commsType == Pstream::nonBlocking)
    {
        label nOutstanding = Pstream::nRequests();

        if (!contiguous<T>())
        {
            PstreamBuffers pBufs(Pstream::nonBlocking, tag);

            // Stream data into buffer
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    // Put data into send buffer
                    UOPstream toDomain(domain, pBufs);
                    toDomain << UIndirectList<T>(field, map);
                }
            }

            // Start receiving. Do not block.
            pBufs.finishedSends(false);

            {
                // Set up 'send' to myself
                List<T> mySubField(field, subMap[Pstream::myProcNo()]);
                // Combine bits. Note that can reuse field storage
                field.setSize(constructSize);
                field = nullValue;
                // Receive sub field from myself
                {
                    const labelList& map = constructMap[Pstream::myProcNo()];

                    forAll(map, i)
                    {
                        cop(field[map[i]], mySubField[i]);
                    }
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

                    forAll(map, i)
                    {
                        cop(field[map[i]], recvField[i]);
                    }
                }
            }
        }
        else
        {
            // Set up sends to neighbours

            List<List<T > > sendFields(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    List<T>& subField = sendFields[domain];
                    subField.setSize(map.size());
                    forAll(map, i)
                    {
                        subField[i] = field[map[i]];
                    }

                    OPstream::write
                    (
                        Pstream::nonBlocking,
                        domain,
                        reinterpret_cast<const char*>(subField.begin()),
                        subField.size()*sizeof(T),
                        tag
                    );
                }
            }

            // Set up receives from neighbours

            List<List<T > > recvFields(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    recvFields[domain].setSize(map.size());
                    UIPstream::read
                    (
                        Pstream::nonBlocking,
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
                    subField[i] = field[map[i]];
                }
            }


            // Combine bits. Note that can reuse field storage

            field.setSize(constructSize);
            field = nullValue;

            // Receive sub field from myself (subField)
            {
                const labelList& map = constructMap[Pstream::myProcNo()];
                const List<T>& subField = sendFields[Pstream::myProcNo()];

                forAll(map, i)
                {
                    cop(field[map[i]], subField[i]);
                }
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

                    forAll(map, i)
                    {
                        cop(field[map[i]], subField[i]);
                    }
                }
            }
        }
    }
    else
    {
        FatalErrorIn("mapDistribute::distribute(..)")
            << "Unknown communication schedule " << commsType
            << abort(FatalError);
    }
}


template<class T>
void Foam::mapDistribute::send(PstreamBuffers& pBufs, const List<T>& field)
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
            toDomain << UIndirectList<T>(field, map);
        }
    }

    // Start sending and receiving but do not block.
    pBufs.finishedSends(false);
}


template<class T>
void Foam::mapDistribute::receive(PstreamBuffers& pBufs, List<T>& field) const
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
                FatalErrorIn
                (
                    "template<class T>\n"
                    "void mapDistribute::receive\n"
                    "(\n"
                    "    PstreamBuffers&,\n"
                    "    List<T>&\n"
                    ")\n"
                )   << "Expected from processor " << domain
                    << " " << map.size() << " but received "
                    << recvField.size() << " elements."
                    << abort(FatalError);
            }

            forAll(map, i)
            {
                field[map[i]] = recvField[i];
            }
        }
    }
}


// In case of no transform: copy elements
template<class T>
void Foam::mapDistribute::applyDummyTransforms(List<T>& field) const
{
    forAll(transformElements_, trafoI)
    {
        const labelList& elems = transformElements_[trafoI];

        label n = transformStart_[trafoI];

        forAll(elems, i)
        {
            field[n++] = field[elems[i]];
        }
    }
}


// In case of no transform: copy elements
template<class T>
void Foam::mapDistribute::applyDummyInverseTransforms(List<T>& field) const
{
    forAll(transformElements_, trafoI)
    {
        const labelList& elems = transformElements_[trafoI];
        label n = transformStart_[trafoI];

        forAll(elems, i)
        {
            field[elems[i]] = field[n++];
        }
    }
}


// Calculate transformed elements.
template<class T, class TransformOp>   //, class CombineOp>
void Foam::mapDistribute::applyTransforms
(
    const globalIndexAndTransform& globalTransforms,
    List<T>& field,
    const TransformOp& top
) const
{
    const List<vectorTensorTransform>& totalTransform =
        globalTransforms.transformPermutations();

    forAll(totalTransform, trafoI)
    {
        const vectorTensorTransform& vt = totalTransform[trafoI];
        const labelList& elems = transformElements_[trafoI];
        label n = transformStart_[trafoI];

        // Could be optimised to avoid memory allocations
        List<T> transformFld(UIndirectList<T>(field, elems));
        top(vt, true, transformFld);

        forAll(transformFld, i)
        {
            //cop(field[n++], transformFld[i]);
            field[n++] = transformFld[i];
        }
    }
}


// Calculate transformed elements.
template<class T, class TransformOp>   //, class CombineOp>
void Foam::mapDistribute::applyInverseTransforms
(
    const globalIndexAndTransform& globalTransforms,
    List<T>& field,
    const TransformOp& top
) const
{
    const List<vectorTensorTransform>& totalTransform =
        globalTransforms.transformPermutations();

    forAll(totalTransform, trafoI)
    {
        const vectorTensorTransform& vt = totalTransform[trafoI];
        const labelList& elems = transformElements_[trafoI];
        label n = transformStart_[trafoI];

        // Could be optimised to avoid memory allocations
        List<T> transformFld(SubList<T>(field, elems.size(), n));
        top(vt, false, transformFld);

        forAll(transformFld, i)
        {
            //cop(field[elems[i]], transformFld[i]);
            field[elems[i]] = transformFld[i];
        }
    }
}


//- Distribute data using default commsType.
template<class T>
void Foam::mapDistribute::distribute
(
    DynamicList<T>& fld,
    const bool dummyTransform,
    const int tag
) const
{
    fld.shrink();

    List<T>& fldList = static_cast<List<T>& >(fld);

    distribute(fldList, dummyTransform, tag);

    fld.setCapacity(fldList.size());
}


//- Distribute data using default commsType.
template<class T>
void Foam::mapDistribute::distribute
(
    List<T>& fld,
    const bool dummyTransform,
    const int tag
) const
{
    if (Pstream::defaultCommsType == Pstream::nonBlocking)
    {
        distribute
        (
            Pstream::nonBlocking,
            List<labelPair>(),
            constructSize_,
            subMap_,
            constructMap_,
            fld,
            tag
        );
    }
    else if (Pstream::defaultCommsType == Pstream::scheduled)
    {
        distribute
        (
            Pstream::scheduled,
            schedule(),
            constructSize_,
            subMap_,
            constructMap_,
            fld,
            tag
        );
    }
    else
    {
        distribute
        (
            Pstream::blocking,
            List<labelPair>(),
            constructSize_,
            subMap_,
            constructMap_,
            fld,
            tag
        );
    }

    //- Fill in transformed slots with copies
    if (dummyTransform)
    {
        applyDummyTransforms(fld);
    }
}


//- Reverse distribute data using default commsType.
template<class T>
void Foam::mapDistribute::reverseDistribute
(
    const label constructSize,
    List<T>& fld,
    const bool dummyTransform,
    const int tag
) const
{
    if (dummyTransform)
    {
        applyDummyInverseTransforms(fld);
    }

    if (Pstream::defaultCommsType == Pstream::nonBlocking)
    {
        distribute
        (
            Pstream::nonBlocking,
            List<labelPair>(),
            constructSize,
            constructMap_,
            subMap_,
            fld,
            tag
        );
    }
    else if (Pstream::defaultCommsType == Pstream::scheduled)
    {
        distribute
        (
            Pstream::scheduled,
            schedule(),
            constructSize,
            constructMap_,
            subMap_,
            fld,
            tag
        );
    }
    else
    {
        distribute
        (
            Pstream::blocking,
            List<labelPair>(),
            constructSize,
            constructMap_,
            subMap_,
            fld,
            tag
        );
    }
}


//- Reverse distribute data using default commsType.
//  Since constructSize might be larger than supplied size supply
//  a nullValue
template<class T>
void Foam::mapDistribute::reverseDistribute
(
    const label constructSize,
    const T& nullValue,
    List<T>& fld,
    const bool dummyTransform,
    const int tag
) const
{
    if (dummyTransform)
    {
        applyDummyInverseTransforms(fld);
    }

    if (Pstream::defaultCommsType == Pstream::nonBlocking)
    {
        distribute
        (
            Pstream::nonBlocking,
            List<labelPair>(),
            constructSize,
            constructMap_,
            subMap_,
            fld,
            eqOp<T>(),
            nullValue,
            tag
        );
    }
    else if (Pstream::defaultCommsType == Pstream::scheduled)
    {
        distribute
        (
            Pstream::scheduled,
            schedule(),
            constructSize,
            constructMap_,
            subMap_,
            fld,
            eqOp<T>(),
            nullValue,
            tag
        );
    }
    else
    {
        distribute
        (
            Pstream::blocking,
            List<labelPair>(),
            constructSize,
            constructMap_,
            subMap_,
            fld,
            eqOp<T>(),
            nullValue,
            tag
        );
    }
}


//- Distribute data using default commsType.
template<class T, class TransformOp>
void Foam::mapDistribute::distribute
(
    const globalIndexAndTransform& git,
    List<T>& fld,
    const TransformOp& top,
    const int tag
) const
{
    // Distribute. Leave out dummy transforms since we're doing them ourselves
    distribute(fld, false, tag);
    // Do transforms
    applyTransforms(git, fld, top);
}


template<class T, class TransformOp>
void Foam::mapDistribute::reverseDistribute
(
    const globalIndexAndTransform& git,
    const label constructSize,
    List<T>& fld,
    const TransformOp& top,
    const int tag
) const
{
    // Fill slots with reverse-transformed data. Note that it also copies
    // back into the non-remote part of fld even though these values are not
    // used.
    applyInverseTransforms(git, fld, top);

    // And send back (the remote slots). Disable dummy transformations.
    reverseDistribute(constructSize, fld, false, tag);
}


template<class T, class TransformOp>
void Foam::mapDistribute::reverseDistribute
(
    const globalIndexAndTransform& git,
    const label constructSize,
    const T& nullValue,
    List<T>& fld,
    const TransformOp& top,
    const int tag
) const
{
    // Fill slots with reverse-transformed data Note that it also copies
    // back into the non-remote part of fld even though these values are not
    // used.
    applyInverseTransforms(git, fld, top);   //, eqOp<T>());

    // And send back (the remote slots) Disable dummy transformations.
    reverseDistribute(constructSize, nullValue, fld, false, tag);
}


// ************************************************************************* //
