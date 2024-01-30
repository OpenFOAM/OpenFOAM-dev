/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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

#include "nonConformalMappedFvPatchBase.H"
#include "nonConformalMappedPolyFacesFvsPatchLabelField.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::nonConformalMappedFvPatchBase::map
(
    const fvsPatchLabelField& srcPolyFacesPf,
    const Field<Type>& srcFld,
    const fvsPatchLabelField& tgtPolyFacesPf
)
{
    if (!Pstream::parRun())
    {
        // Non-conformal mapped patches are ordered. So just return the field.
        return srcFld;
    }
    else
    {
        // Otherwise this is a little more involved. Use the proc offsets and
        // sizes stored in the polyFacesPf fields to send and receive in blocks.

        typedef
            nonConformalMappedPolyFacesFvsPatchLabelField
            ncmpfFvsPatchLabelField;
        const ncmpfFvsPatchLabelField& srcNcmPolyFacesPf =
            refCast<const ncmpfFvsPatchLabelField>(srcPolyFacesPf);
        const ncmpfFvsPatchLabelField& tgtNcmPolyFacesPf =
            refCast<const ncmpfFvsPatchLabelField>(tgtPolyFacesPf);
        const labelList& srcProcOffsets = srcNcmPolyFacesPf.procOffsets();
        const labelList srcProcSizes(srcNcmPolyFacesPf.procSizes());
        const labelList& tgtProcOffsets = tgtNcmPolyFacesPf.procOffsets();
        const labelList tgtProcSizes(tgtNcmPolyFacesPf.procSizes());

        tmp<Field<Type>> tTgtFld
        (
            new Field<Type>(tgtProcOffsets.last() + tgtProcSizes.last())
        );
        Field<Type>& tgtFld = tTgtFld.ref();

        // Do local copy
        {
            SubList<Type> srcProcFld
            (
                srcFld,
                srcProcSizes[Pstream::myProcNo()],
                srcProcOffsets[Pstream::myProcNo()]
            );

            SubList<Type> tgtProcFld
            (
                tgtFld,
                tgtProcSizes[Pstream::myProcNo()],
                tgtProcOffsets[Pstream::myProcNo()]
            );

            tgtProcFld = srcProcFld;
        }

        // Do remote transfers
        {
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

            forAll(srcProcOffsets, proci)
            {
                if (proci != Pstream::myProcNo())
                {
                    SubList<Type> srcProcFld
                    (
                        srcFld,
                        srcProcSizes[proci],
                        srcProcOffsets[proci]
                    );

                    UOPstream(proci, pBufs)() << srcProcFld;
                }
            }

            pBufs.finishedSends();

            forAll(tgtProcOffsets, proci)
            {
                if (proci != Pstream::myProcNo())
                {
                    SubList<Type> tgtProcFld
                    (
                        tgtFld,
                        tgtProcSizes[proci],
                        tgtProcOffsets[proci]
                    );

                    UIPstream(proci, pBufs)() >> tgtProcFld;
                }
            }
        }

        return tTgtFld;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::nonConformalMappedFvPatchBase::map
(
    const fvsPatchLabelField& srcPolyFacesPf,
    Field<Type>&& srcFld,
    const fvsPatchLabelField& tgtPolyFacesPf
)
{
    if (!Pstream::parRun())
    {
        // Non-conformal mapped patches are ordered. So just return the field.
        return tmp<Field<Type>>(new Field<Type>(std::move(srcFld)));
    }
    else
    {
        // Defer to the above to communicate and construct the result
        tmp<Field<Type>> tResult =
            map(srcPolyFacesPf, srcFld, tgtPolyFacesPf);
        return tResult;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::nonConformalMappedFvPatchBase::map
(
    const fvsPatchLabelField& srcPolyFacesPf,
    const tmp<Field<Type>>& srcFld,
    const fvsPatchLabelField& tgtPolyFacesPf
)
{
    if (!Pstream::parRun())
    {
        // Non-conformal mapped patches are ordered. So just return the field.
        return srcFld;
    }
    else
    {
        // Defer to the above to communicate and construct the result
        tmp<Field<Type>> tResult =
            map(srcPolyFacesPf, srcFld(), tgtPolyFacesPf);
        srcFld.clear();
        return tResult;
    }
}


// ************************************************************************* //
