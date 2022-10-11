/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2022 OpenFOAM Foundation
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

#include "meshToMesh.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::meshToMesh::add
(
    UList<Type>& fld,
    const label offset
)
{
    forAll(fld, i)
    {
        fld[i] += offset;
    }
}


template<class Type>
void Foam::meshToMesh::mapSrcToTgt
(
    const UList<Type>& srcField,
    List<Type>& result
) const
{
    if (result.size() != tgtLocalSrcCells_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to target mesh size" << nl
            << "    source mesh    = " << srcLocalTgtCells_.size() << nl
            << "    target mesh    = " << tgtLocalSrcCells_.size() << nl
            << "    supplied field = " << result.size()
            << abort(FatalError);
    }

    if (!isSingleProcess())
    {
        const distributionMap& map = srcMapPtr_();

        List<Type> work(srcField);
        map.distribute(work);

        forAll(result, tgtCelli)
        {
            if (tgtLocalSrcCells_[tgtCelli].size())
            {
                result[tgtCelli] *= (1.0 - sum(tgtWeights_[tgtCelli]));
                forAll(tgtLocalSrcCells_[tgtCelli], i)
                {
                    result[tgtCelli] +=
                        tgtWeights_[tgtCelli][i]
                       *work[tgtLocalSrcCells_[tgtCelli][i]];
                }
            }
        }
    }
    else
    {
        forAll(result, tgtCelli)
        {
            if (tgtLocalSrcCells_[tgtCelli].size())
            {
                result[tgtCelli] *= (1.0 - sum(tgtWeights_[tgtCelli]));
                forAll(tgtLocalSrcCells_[tgtCelli], i)
                {
                    result[tgtCelli] +=
                        tgtWeights_[tgtCelli][i]
                       *srcField[tgtLocalSrcCells_[tgtCelli][i]];
                }
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::meshToMesh::mapSrcToTgt
(
    const Field<Type>& srcField
) const
{
    tmp<Field<Type>> tresult
    (
        new Field<Type>
        (
            tgtLocalSrcCells_.size(),
            Zero
        )
    );

    mapSrcToTgt(srcField, tresult.ref());

    return tresult;
}


template<class Type>
void Foam::meshToMesh::mapTgtToSrc
(
    const UList<Type>& tgtField,
    List<Type>& result
) const
{
    if (result.size() != srcLocalTgtCells_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to source mesh size" << nl
            << "    source mesh    = " << srcLocalTgtCells_.size() << nl
            << "    target mesh    = " << tgtLocalSrcCells_.size() << nl
            << "    supplied field = " << result.size()
            << abort(FatalError);
    }

    if (!isSingleProcess())
    {
        const distributionMap& map = tgtMapPtr_();

        List<Type> work(tgtField);
        map.distribute(work);

        forAll(result, srcCelli)
        {
            if (srcLocalTgtCells_[srcCelli].size())
            {
                result[srcCelli] *= (1.0 - sum(srcWeights_[srcCelli]));
                forAll(srcLocalTgtCells_[srcCelli], i)
                {
                    result[srcCelli] +=
                        srcWeights_[srcCelli][i]
                       *work[srcLocalTgtCells_[srcCelli][i]];
                }
            }
        }
    }
    else
    {
        forAll(result, srcCelli)
        {
            if (srcLocalTgtCells_[srcCelli].size())
            {
                result[srcCelli] *= (1.0 - sum(srcWeights_[srcCelli]));
                forAll(srcLocalTgtCells_[srcCelli], i)
                {
                    result[srcCelli] +=
                        srcWeights_[srcCelli][i]
                       *tgtField[srcLocalTgtCells_[srcCelli][i]];
                }
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::meshToMesh::mapTgtToSrc
(
    const Field<Type>& tgtField
) const
{
    tmp<Field<Type>> tresult
    (
        new Field<Type>
        (
            srcLocalTgtCells_.size(),
            Zero
        )
    );

    mapTgtToSrc(tgtField, tresult.ref());

    return tresult;
}


// ************************************************************************* //
