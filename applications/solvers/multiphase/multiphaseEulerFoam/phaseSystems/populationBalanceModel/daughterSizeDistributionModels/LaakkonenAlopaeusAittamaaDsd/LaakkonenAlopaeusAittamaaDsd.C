/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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

#include "LaakkonenAlopaeusAittamaaDsd.H"
#include "addToRunTimeSelectionTable.H"
#include "breakupModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace daughterSizeDistributionModels
{
    defineTypeNameAndDebug(LaakkonenAlopaeusAittamaaDsd, 0);
    addToRunTimeSelectionTable
    (
        daughterSizeDistributionModel,
        LaakkonenAlopaeusAittamaaDsd,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::daughterSizeDistributionModels::
LaakkonenAlopaeusAittamaaDsd::LaakkonenAlopaeusAittamaaDsd
(
    const breakupModel& breakup,
    const dictionary& dict
)
:
    daughterSizeDistributionModel(breakup, dict),
    C4_(dimensionedScalar::lookupOrDefault("C4", dict, dimless, 18.25))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::daughterSizeDistributionModels::
LaakkonenAlopaeusAittamaaDsd::~LaakkonenAlopaeusAittamaaDsd()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar
Foam::diameterModels::daughterSizeDistributionModels::
LaakkonenAlopaeusAittamaaDsd::antiderivative
(
    const dimensionedScalar& xk,
    const dimensionedScalar& v,
    const dimensionedScalar& bndr,
    const dimensionedScalar range
) const
{
    return
        (4.0/3.0 + C4_/3)
       *(
            pow(xk, -C4_ - 3)*pow(xk - v, C4_)*(v - xk)
           *(
                (C4_ + 1)*(C4_ + 2)*(C4_ + 3)*pow3(v)
              - (C4_ + 1)*(C4_ + 2)*(bndr*(C4_ + 4) - 3*xk)*sqr(v)
              - 2*v*xk*(C4_ + 1)*(bndr*(C4_ + 4) - 3*xk)
              - 2*bndr*C4_*sqr(xk)
              + 6*pow3(xk)
              - 8*bndr*sqr(xk)
            )
        )/(2*range*(C4_ + 4));
}


Foam::dimensionedScalar
Foam::diameterModels::daughterSizeDistributionModels::
LaakkonenAlopaeusAittamaaDsd::calcNik
(
    const label i,
    const label k
) const
{
    const dimensionedScalar& x0 = breakup_.popBal().sizeGroups()[0].x();
    dimensionedScalar xi = breakup_.popBal().sizeGroups()[i].x() - x0;
    dimensionedScalar xk = breakup_.popBal().sizeGroups()[k].x() - x0;
    const UPtrList<sizeGroup>& sizeGroups = breakup_.popBal().sizeGroups();

    if (i == 0)
    {
        dimensionedScalar xii = sizeGroups[i+1].x() - x0;

        if (k == 0)
        {
            return 1.0;
        }

        return
            antiderivative(xk, xi, xii, (xii-xi))
          - antiderivative(xk, xii, xii, (xii-xi));
    }
    else if (i == k)
    {
        dimensionedScalar x = sizeGroups[i-1].x() - x0;

        return
            antiderivative(xk, xi, x, (xi-x))
          - antiderivative(xk, x, x, (xi-x));
    }
    else
    {
        dimensionedScalar x = sizeGroups[i-1].x() - x0;
        dimensionedScalar xii = sizeGroups[i+1].x() - x0;

        return
            antiderivative(xk, xi, xii, (xii-xi))
          - antiderivative(xk, xii, xii, (xii-xi))
          + antiderivative(xk, xi, x, (xi-x))
          - antiderivative(xk, x, x, (xi-x));
    }
}


// ************************************************************************* //
