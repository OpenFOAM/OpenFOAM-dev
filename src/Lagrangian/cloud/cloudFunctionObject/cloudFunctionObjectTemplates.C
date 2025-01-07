/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "cloudFunctionObject.H"
#include "cloud.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cloud>
bool Foam::functionObjects::cloudFunctionObject::isCloud() const
{
    return isA<Cloud>(cloud_);
}


template<class Cloud>
const Cloud& Foam::functionObjects::cloudFunctionObject::cloud() const
{
    if (!isA<Cloud>(cloud_))
    {
        FatalErrorInFunction
            << "The cloud function object '" << function_.name()
            << "' of cloud '" << cloud_.mesh().name()
            << "' requires a cloud of type '" << Cloud::typeName
            << "' (or a derivation thereof), rather than '"
            << cloud_.type() << "'" << exit(FatalError);
    }

    return refCast<const Cloud>(cloud_);
}


// ************************************************************************* //
