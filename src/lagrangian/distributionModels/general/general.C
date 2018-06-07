/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "general.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
    defineTypeNameAndDebug(general, 0);
    addToRunTimeSelectionTable(distributionModel, general, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::general::general
(
    const dictionary& dict,
    Random& rndGen
)
:
    distributionModel(typeName, dict, rndGen),
    xy_(distributionModelDict_.lookup("distribution")),
    nEntries_(xy_.size()),
    minValue_(xy_[0][0]),
    maxValue_(xy_[nEntries_-1][0]),
    meanValue_(0.0),
    integral_(nEntries_)
{
    check();

    // normalize the cumulative distributionModel

    integral_[0] = 0.0;
    for (label i=1; i<nEntries_; i++)
    {

        scalar k = (xy_[i][1] - xy_[i-1][1])/(xy_[i][0] - xy_[i-1][0]);
        scalar d = xy_[i-1][1] - k*xy_[i-1][0];
        scalar y1 = xy_[i][0]*(0.5*k*xy_[i][0] + d);
        scalar y0 = xy_[i-1][0]*(0.5*k*xy_[i-1][0] + d);
        scalar area = y1 - y0;

        integral_[i] = area + integral_[i-1];
    }

    scalar sumArea = integral_.last();

    meanValue_ = sumArea/(maxValue_ - minValue_);

    for (label i=0; i<nEntries_; i++)
    {
        xy_[i][1] /= sumArea;
        integral_[i] /= sumArea;
    }

}


Foam::distributionModels::general::general(const general& p)
:
    distributionModel(p),
    xy_(p.xy_),
    nEntries_(p.nEntries_),
    minValue_(p.minValue_),
    maxValue_(p.maxValue_),
    integral_(p.integral_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributionModels::general::~general()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::general::sample() const
{
    scalar y = rndGen_.sample01<scalar>();

    // find the interval where y is in the table
    label n=1;
    while (integral_[n] <= y)
    {
        n++;
    }

    scalar k = (xy_[n][1] - xy_[n-1][1])/(xy_[n][0] - xy_[n-1][0]);
    scalar d = xy_[n-1][1] - k*xy_[n-1][0];

    scalar alpha = y + xy_[n-1][0]*(0.5*k*xy_[n-1][0] + d) - integral_[n-1];
    scalar x = 0.0;

    // if k is small it is a linear equation, otherwise it is of second order
    if (mag(k) > small)
    {
        scalar p = 2.0*d/k;
        scalar q = -2.0*alpha/k;
        scalar sqrtEr = sqrt(0.25*p*p - q);

        scalar x1 = -0.5*p + sqrtEr;
        scalar x2 = -0.5*p - sqrtEr;
        if ((x1 >= xy_[n-1][0]) && (x1 <= xy_[n][0]))
        {
            x = x1;
        }
        else
        {
            x = x2;
        }
    }
    else
    {
        x = alpha/d;
    }

    return x;
}


Foam::scalar Foam::distributionModels::general::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::distributionModels::general::maxValue() const
{
    return maxValue_;
}


Foam::scalar Foam::distributionModels::general::meanValue() const
{
    return meanValue_;
}


// ************************************************************************* //
