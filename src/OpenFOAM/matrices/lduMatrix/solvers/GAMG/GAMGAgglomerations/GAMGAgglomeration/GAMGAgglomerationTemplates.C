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

\*---------------------------------------------------------------------------*/

#include "GAMGAgglomeration.H"
#include "mapDistribute.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::GAMGAgglomeration::gatherList
(
    const label comm,
    const labelList& procIDs,
    const Type& myVal,
    List<Type>& allVals,
    const int tag
)
{
    if (Pstream::myProcNo(comm) == procIDs[0])
    {
        allVals.setSize(procIDs.size());

        allVals[0] = myVal;
        for (label i=1; i<procIDs.size(); i++)
        {
            IPstream fromSlave
            (
                Pstream::commsTypes::scheduled,
                procIDs[i],
                0,
                tag,
                comm
            );

            fromSlave >> allVals[i];
        }
    }
    else
    {
        OPstream toMaster
        (
            Pstream::commsTypes::scheduled,
            procIDs[0],
            0,
            tag,
            comm
        );
        toMaster << myVal;
    }
}


template<class Type>
void Foam::GAMGAgglomeration::restrictField
(
    Field<Type>& cf,
    const Field<Type>& ff,
    const labelList& fineToCoarse
) const
{
    cf = Zero;

    forAll(ff, i)
    {
        cf[fineToCoarse[i]] += ff[i];
    }
}


template<class Type>
void Foam::GAMGAgglomeration::restrictField
(
    Field<Type>& cf,
    const Field<Type>& ff,
    const label fineLevelIndex,
    const bool procAgglom
) const
{
    const labelList& fineToCoarse = restrictAddressing_[fineLevelIndex];

    if (!procAgglom && ff.size() != fineToCoarse.size())
    {
        FatalErrorInFunction
            << "field does not correspond to level " << fineLevelIndex
            << " sizes: field = " << ff.size()
            << " level = " << fineToCoarse.size()
            << abort(FatalError);
    }

    restrictField(cf, ff, fineToCoarse);

    label coarseLevelIndex = fineLevelIndex+1;

    if (procAgglom && hasProcMesh(coarseLevelIndex))
    {
        label fineComm = UPstream::parent(procCommunicator_[coarseLevelIndex]);

        const List<label>& procIDs = agglomProcIDs(coarseLevelIndex);
        const labelList& offsets = cellOffsets(coarseLevelIndex);

        globalIndex::gather
        (
            offsets,
            fineComm,
            procIDs,
            cf,
            UPstream::msgType(),
            Pstream::commsTypes::nonBlocking // Pstream::commsTypes::scheduled
        );
    }
}


template<class Type>
void Foam::GAMGAgglomeration::restrictFaceField
(
    Field<Type>& cf,
    const Field<Type>& ff,
    const label fineLevelIndex
) const
{
    const labelList& fineToCoarse = faceRestrictAddressing_[fineLevelIndex];

    if (ff.size() != fineToCoarse.size())
    {
        FatalErrorInFunction
            << "field does not correspond to level " << fineLevelIndex
            << " sizes: field = " << ff.size()
            << " level = " << fineToCoarse.size()
            << abort(FatalError);
    }

    cf = Zero;

    forAll(fineToCoarse, ffacei)
    {
        label cFace = fineToCoarse[ffacei];

        if (cFace >= 0)
        {
            cf[cFace] += ff[ffacei];
        }
    }
}


template<class Type>
void Foam::GAMGAgglomeration::prolongField
(
    Field<Type>& ff,
    const Field<Type>& cf,
    const label levelIndex,
    const bool procAgglom
) const
{
    const labelList& fineToCoarse = restrictAddressing_[levelIndex];

    label coarseLevelIndex = levelIndex+1;

    if (procAgglom && hasProcMesh(coarseLevelIndex))
    {
        label coarseComm = UPstream::parent
        (
            procCommunicator_[coarseLevelIndex]
        );

        const List<label>& procIDs = agglomProcIDs(coarseLevelIndex);
        const labelList& offsets = cellOffsets(coarseLevelIndex);

        label localSize = nCells_[levelIndex];

        Field<Type> allCf(localSize);
        globalIndex::scatter
        (
            offsets,
            coarseComm,
            procIDs,
            cf,
            allCf,
            UPstream::msgType(),
            Pstream::commsTypes::nonBlocking // Pstream::commsTypes::scheduled
        );

        forAll(fineToCoarse, i)
        {
            ff[i] = allCf[fineToCoarse[i]];
        }
    }
    else
    {
        forAll(fineToCoarse, i)
        {
            ff[i] = cf[fineToCoarse[i]];
        }
    }
}


// ************************************************************************* //
