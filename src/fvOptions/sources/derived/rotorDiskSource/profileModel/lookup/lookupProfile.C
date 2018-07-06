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

#include "lookupProfile.H"
#include "addToRunTimeSelectionTable.H"
#include "vector.H"
#include "unitConversion.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lookupProfile, 0);
    addToRunTimeSelectionTable(profileModel, lookupProfile, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::lookupProfile::interpolateWeights
(
    const scalar& xIn,
    const List<scalar>& values,
    label& i1,
    label& i2,
    scalar& ddx
) const
{
    i2 = 0;
    label nElem = values.size();

    if (nElem == 1)
    {
        i1 = i2;
        ddx = 0.0;
        return;
    }
    else
    {
        while ((i2 < nElem) && (values[i2] < xIn))
        {
            i2++;
        }

        if (i2 == 0)
        {
            i1 = i2;
            ddx = 0.0;
            return;
        }
        else if (i2 == nElem)
        {
            i2 = nElem - 1;
            i1 = i2;
            ddx = 0.0;
            return;
        }
        else
        {
            i1 = i2 - 1;
            ddx = (xIn - values[i1])/(values[i2] - values[i1]);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lookupProfile::lookupProfile
(
    const dictionary& dict,
    const word& modelName
)
:
    profileModel(dict, modelName),
    AOA_(),
    Cd_(),
    Cl_()
{
    List<vector> data;
    if (readFromFile())
    {
        IFstream is(fName_);
        is  >> data;
    }
    else
    {
        dict.lookup("data") >> data;
    }

    if (data.size() > 0)
    {
        AOA_.setSize(data.size());
        Cd_.setSize(data.size());
        Cl_.setSize(data.size());

        forAll(data, i)
        {
            AOA_[i] = degToRad(data[i][0]);
            Cd_[i] = data[i][1];
            Cl_[i] = data[i][2];
        }
    }
    else
    {
        FatalErrorInFunction
            << "No profile data specified" << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lookupProfile::Cdl(const scalar alpha, scalar& Cd, scalar& Cl) const
{
    label i1 = -1;
    label i2 = -1;
    scalar invAlpha = -1.0;
    interpolateWeights(alpha, AOA_, i1, i2, invAlpha);

    Cd = invAlpha*(Cd_[i2] - Cd_[i1]) + Cd_[i1];
    Cl = invAlpha*(Cl_[i2] - Cl_[i1]) + Cl_[i1];
}


// ************************************************************************* //
