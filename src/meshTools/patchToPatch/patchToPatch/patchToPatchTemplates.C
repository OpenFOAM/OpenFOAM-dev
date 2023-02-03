/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "patchToPatch.H"
#include "patchToPatchTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchToPatch::srcToTgt(const Field<Type>& srcFld) const
{
    return
        patchToPatchTools::interpolate
        (
            tgtLocalSrcFaces_,
            tgtWeights()(),
            srcMapPtr_,
            srcFld
        );
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::patchToPatch::srcToTgt
(
    const Field<Type>& srcFld,
    const Field<Type>& leftOverTgtFld
) const
{
    return
        patchToPatchTools::interpolate
        (
            tgtLocalSrcFaces_,
            tgtWeights()(),
            srcMapPtr_,
            srcFld,
            leftOverTgtFld
        );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchToPatch::tgtToSrc(const Field<Type>& tgtFld) const
{
    return
        patchToPatchTools::interpolate
        (
            srcLocalTgtFaces_,
            srcWeights()(),
            tgtMapPtr_,
            tgtFld
        );
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::patchToPatch::tgtToSrc
(
    const Field<Type>& tgtFld,
    const Field<Type>& leftOverSrcFld
) const
{
    return
        patchToPatchTools::interpolate
        (
            srcLocalTgtFaces_,
            srcWeights()(),
            tgtMapPtr_,
            tgtFld,
            leftOverSrcFld
        );
}


// ************************************************************************* //
