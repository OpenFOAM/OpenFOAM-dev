/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

#include "uniformBinaryDsd.H"
#include "addToRunTimeSelectionTable.H"
#include "breakupModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace daughterSizeDistributionModels
{
    defineTypeNameAndDebug(uniformBinaryDsd, 0);
    addToRunTimeSelectionTable
    (
        daughterSizeDistributionModel,
        uniformBinaryDsd,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::daughterSizeDistributionModels::uniformBinaryDsd::
uniformBinaryDsd
(
    const breakupModel& breakup,
    const dictionary& dict
)
:
    daughterSizeDistributionModel(breakup, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::daughterSizeDistributionModels::uniformBinaryDsd::
~uniformBinaryDsd()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar
Foam::diameterModels::daughterSizeDistributionModels::uniformBinaryDsd::n
(
    const label i,
    const label k
) const
{
    const sizeGroup& fi = *breakup_.popBal().sizeGroups()[i];
    const sizeGroup& fk = *breakup_.popBal().sizeGroups()[k];

    const dimensionedScalar& xi = fi.x();
    const dimensionedScalar& xk = fk.x();
    dimensionedScalar x("x", dimVolume, Zero);
    dimensionedScalar xii("xii", dimVolume, Zero);

    if (i == 0)
    {
       x = xi;
    }
    else
    {
       x = breakup_.popBal().sizeGroups()[i-1]->x();
    }

    if (i == breakup_.popBal().sizeGroups().size() - 1)
    {
        xii = fi.x();
    }
    else
    {
        xii = breakup_.popBal().sizeGroups()[i+1]->x();
    }

    if (i == 0)
    {
        return
            1.0/(xii - xi)
           *(
                (xii*xii - 0.5*sqr(xii))
              - (xi*xii - 0.5*sqr(xi))
            )
           *2.0/xk;
    }
    else if (i == k)
    {
        return
            1.0/(xi - x)
           *(
                (0.5*sqr(xi) - xi*x)
              - (0.5*sqr(x) - x*x)
            )
           *2.0/xk;
    }
    else
    {
        return
            1.0/(xii - xi)
           *(
                (xii*xii - 0.5*sqr(xii))
              - (xi*xii - 0.5*sqr(xi))
            )
           *2.0/xk
          + 1.0/(xi - x)
           *(
                (0.5*sqr(xi) - xi*x)
              - (0.5*sqr(x) - x*x)
            )
           *2.0/xk;
    }
}


// ************************************************************************* //
