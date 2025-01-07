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

#include "cloudLagrangianModel.H"
#include "cloud.H"
#include "CloudTypes.H"

// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

template<class Cloud, class ... Clouds>
bool Foam::cloudLagrangianModel::isCloud() const
{
    return CloudTypes<Cloud, Clouds ...>::isA(cloud_);
}


template<class Cloud, class ... Clouds>
void Foam::cloudLagrangianModel::assertCloud() const
{
    if (!isCloud<Cloud, Clouds ...>())
    {
        FatalErrorInFunction
            << "The Larangian model '" << name_ << "' of cloud '"
            << cloud_.mesh().name() << "' requires a cloud of type "
            << CloudTypes<Cloud, Clouds ...>::typesString("or").c_str()
            << " (or a derivation thereof), rather than '" << cloud_.type()
            << "'" << exit(FatalError);
    }
}


template<class Cloud>
const Cloud& Foam::cloudLagrangianModel::cloud() const
{
    assertCloud<Cloud>();
    return refCast<const Cloud>(cloud_);
}


// ************************************************************************* //
