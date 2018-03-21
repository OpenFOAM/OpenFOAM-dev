/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
#include "IFstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::fileOperations::masterUncollatedFileOperation::scatterList
(
    const UList<Type>& masterLst,
    const int tag,
    const label comm
) const
{
    // TBD: more efficient scatter
    PstreamBuffers pBufs(UPstream::commsTypes::nonBlocking, tag, comm);
    if (Pstream::master(comm))
    {
        for (label proci = 1; proci < Pstream::nProcs(comm); proci++)
        {
            UOPstream os(proci, pBufs);
            os << masterLst[proci];
        }
    }
    pBufs.finishedSends();

    Type myResult;

    if (Pstream::master(comm))
    {
        myResult = masterLst[Pstream::myProcNo(comm)];
    }
    else
    {
        UIPstream is(Pstream::masterNo(), pBufs);
        is >> myResult;
    }
    return myResult;
}


template<class Type, class fileOp>
Type Foam::fileOperations::masterUncollatedFileOperation::masterOp
(
    const fileName& fName,
    const fileOp& fop,
    const int tag,
    const label comm
) const
{
    if (IFstream::debug)
    {
        Pout<< "masterUncollatedFileOperation::masterOp : Operation "
            << typeid(fileOp).name()
            << " on " << fName << endl;
    }
    if (Pstream::parRun())
    {
        List<fileName> filePaths(Pstream::nProcs(comm));
        filePaths[Pstream::myProcNo(comm)] = fName;
        Pstream::gatherList(filePaths, tag, comm);

        List<Type> result(filePaths.size());
        if (Pstream::master(comm))
        {
            result = fop(filePaths[0]);
            for (label i = 1; i < filePaths.size(); i++)
            {
                if (filePaths[i] != filePaths[0])
                {
                    result[i] = fop(filePaths[i]);
                }
            }
        }

        return scatterList(result, tag, comm);
    }
    else
    {
        return fop(fName);
    }
}


template<class Type, class fileOp>
Type Foam::fileOperations::masterUncollatedFileOperation::masterOp
(
    const fileName& src,
    const fileName& dest,
    const fileOp& fop,
    const int tag,
    const label comm
) const
{
    if (IFstream::debug)
    {
        Pout<< "masterUncollatedFileOperation : Operation on src:" << src
            << " dest:" << dest << endl;
    }
    if (Pstream::parRun())
    {
        List<fileName> srcs(Pstream::nProcs(comm));
        srcs[Pstream::myProcNo(comm)] = src;
        Pstream::gatherList(srcs, tag, comm);

        List<fileName> dests(srcs.size());
        dests[Pstream::myProcNo(comm)] = dest;
        Pstream::gatherList(dests, tag, comm);

        List<Type> result(Pstream::nProcs(comm));
        if (Pstream::master(comm))
        {
            result = fop(srcs[0], dests[0]);
            for (label i = 1; i < srcs.size(); i++)
            {
                if (srcs[i] != srcs[0])
                {
                    result[i] = fop(srcs[i], dests[i]);
                }
            }
        }

        return scatterList(result, tag, comm);
    }
    else
    {
        return fop(src, dest);
    }
}


// ************************************************************************* //
