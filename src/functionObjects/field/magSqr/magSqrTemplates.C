/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::magSqr::calcMagSqr()
{
    if (foundObject<VolField<Type>>(fieldName_))
    {
        store
        (
            resultName_,
            Foam::magSqr(lookupObject<VolField<Type>>(fieldName_))
        );

        return true;
    }
    else if (foundObject<SurfaceField<Type>>(fieldName_))
    {
        store
        (
            resultName_,
            Foam::magSqr(lookupObject<SurfaceField<Type>>(fieldName_))
        );

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
