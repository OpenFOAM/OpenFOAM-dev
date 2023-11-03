/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "surfaceToVolVelocity.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

const Foam::volVectorField& Foam::surfaceToVolVelocity
(
    const surfaceVectorField& Uf
)
{
    if (!isFaceVelocity(Uf))
    {
        return volVectorField::null();
    }

    const word UfName(Uf.member());

    // Find where the old-time suffixes end
    string::size_type i = UfName.size();
    while (i > 0 && UfName(i - 2, 2) == "_0") i -= 2;

    // Split into current name and old-time suffix
    const word UfCurName = UfName(i);
    const word Uf_0 = UfName(i, UfName.size());

    if (UfCurName.back() != 'f')
    {
        return volVectorField::null();
    }

    const word Uname =
        IOobject::groupName
        (
            UfCurName(UfCurName.size() - 1) + Uf_0,
            Uf.group()
        );

    if (!Uf.mesh().foundObject<volVectorField>(Uname))
    {
        return volVectorField::null();
    }

    return Uf.mesh().lookupObject<volVectorField>(Uname);
}


// ************************************************************************* //
