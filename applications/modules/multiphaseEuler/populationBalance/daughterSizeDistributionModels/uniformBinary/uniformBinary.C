/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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

#include "uniformBinary.H"
#include "addToRunTimeSelectionTable.H"
#include "breakupModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace daughterSizeDistributionModels
{
    defineTypeNameAndDebug(uniformBinary, 0);
    addToRunTimeSelectionTable
    (
        daughterSizeDistributionModel,
        uniformBinary,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::daughterSizeDistributionModels::uniformBinary::
uniformBinary
(
    const breakupModel& breakup,
    const dictionary& dict
)
:
    daughterSizeDistributionModel(breakup, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalance::daughterSizeDistributionModels::uniformBinary::
~uniformBinary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar
Foam::populationBalance::daughterSizeDistributionModels::uniformBinary::calcNik
(
    const label i,
    const label k
) const
{
    const populationBalanceModel& popBal = breakup_.popBal();

    const dimensionedScalar& x0 = popBal.v(0);
    const dimensionedScalar& xi = popBal.v(i);
    const dimensionedScalar& xk = popBal.v(k);

    if (i == 0)
    {
        if (k == 0)
        {
            return 1;
        }

        return (popBal.v(i+1) - xi)/(xk - x0);
    }
    else if (i == k)
    {
        return (xi - popBal.v(i-1))/(xk - x0);
    }
    else
    {
        return (popBal.v(i+1) - popBal.v(i-1))/(xk - x0);
    }
}


// ************************************************************************* //
