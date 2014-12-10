/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description
    coupledFvPatch is an abstract base class for patches that couple regions
    of the computational domain e.g. cyclic and processor-processor links.

\*---------------------------------------------------------------------------*/

#include "regionCoupledFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionCoupledFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, regionCoupledFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::regionCoupledFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::regionCoupledFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    if (neighbFvPatch().sameRegion())
    {
        return neighbFvPatch().patchInternalField(iF);
    }
    else
    {
        /*
        WarningIn
        (
            "regionCoupledFvPatch::internalFieldTransfer"
            "( const Pstream::commsTypes, const labelUList&)"
            " the internal field can not be transfered "
            " as the neighbFvPatch are in different meshes "
        );
        */
        return tmp<labelField>(new labelField(iF.size(), 0));

    }
    return tmp<labelField>(NULL);
}


// ************************************************************************* //
