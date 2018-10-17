/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    daughterSizeDistributionModel(breakup, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::daughterSizeDistributionModels::
LaakkonenAlopaeusAittamaaDsd::~LaakkonenAlopaeusAittamaaDsd()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar
Foam::diameterModels::daughterSizeDistributionModels::
LaakkonenAlopaeusAittamaaDsd::calcNik
(
    const label i,
    const label k
) const
{
    const dimensionedScalar& xi = breakup_.popBal().sizeGroups()[i].x();
    const dimensionedScalar& xk = breakup_.popBal().sizeGroups()[k].x();
    const UPtrList<sizeGroup>& sizeGroups = breakup_.popBal().sizeGroups();

    if (i == 0)
    {
        const dimensionedScalar& xii = sizeGroups[i+1].x();

        return
            (
                5.0*pow4(xii)*sqr(xk) - 12.0*pow5(xi)*xii + 10.0*pow6(xi)
              - (20.0*pow3(xi)*xii - 15.0*pow4(xi))*sqr(xk) - 6.0*pow5(xii)*xk
              - (24*pow5(xi) - 30*pow4(xi)*xii)*xk + 2*pow6(xii)
            )
           /((xii - xi)*pow5(xk));
    }
    else if (i == k)
    {
        const dimensionedScalar& x = sizeGroups[i-1].x();

        return
            (
                (15.0*pow4(xi) - 20.0*x*pow3(xi))*sqr(xk)
              + 5.0*pow4(x)*sqr(xk) + (30.0*x*pow4(xi) - 24.0*pow5(xi))*xk
              - 6.0*pow5(x)*xk + 10.0*pow6(xi) - 12.0*x*pow5(xi) + 2.0*pow6(x)
            )
           /((xi - x)*pow5(xk));
    }
    else
    {
        const dimensionedScalar& x = sizeGroups[i-1].x();
        const dimensionedScalar& xii = sizeGroups[i+1].x();

        return
            (
                5.0*pow4(xii)*sqr(xk) - 12.0*pow5(xi)*xii + 10.0*pow6(xi)
              - (20.0*pow3(xi)*xii - 15.0*pow4(xi))*sqr(xk) - 6.0*pow5(xii)*xk
              - (24*pow5(xi) - 30*pow4(xi)*xii)*xk + 2*pow6(xii)
            )
           /((xii - xi)*pow5(xk))
          + (
                (15.0*pow4(xi) - 20.0*x*pow3(xi))*sqr(xk)
              + 5.0*pow4(x)*sqr(xk) + (30.0*x*pow4(xi) - 24.0*pow5(xi))*xk
              - 6.0*pow5(x)*xk + 10.0*pow6(xi) - 12.0*x*pow5(xi) + 2.0*pow6(x)
            )
           /((xi - x)*pow5(xk));
    }
}


// ************************************************************************* //
