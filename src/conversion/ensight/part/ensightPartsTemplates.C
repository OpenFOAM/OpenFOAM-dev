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
    Template to write generalized field components

\*---------------------------------------------------------------------------*/

#include "ensightParts.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::ensightParts::writeField
(
    ensightFile& os,
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    // find offset to patch parts (ie, the first face data)
    label patchOffset = 0;
    forAll(partsList_, partI)
    {
        if (partsList_[partI].isFaceData())
        {
            patchOffset = partI;
            break;
        }
    }

    forAll(partsList_, partI)
    {
        label patchI = partI - patchOffset;

        if (partsList_[partI].isCellData())
        {
            partsList_[partI].writeField
            (
                os,
                field
            );
        }
        else if (patchI < field.boundaryField().size())
        {
            partsList_[partI].writeField
            (
                os,
                field.boundaryField()[patchI]
            );
        }
    }
}


// ************************************************************************* //
