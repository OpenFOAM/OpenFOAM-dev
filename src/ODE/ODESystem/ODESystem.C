/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "ODESystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ODESystem::ODESystem()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ODESystem::~ODESystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::ODESystem::check
(
    const scalar x,
    const scalarField& y,
    const scalarField& dy,
    const label li
) const
{
    // Evaluate the derivatives using the derivatives method
    scalarField dfdx0(nEqns());
    derivatives(x, y, li, dfdx0);

    // Evaluate the derivatives and the Jacobian using the Jacobian method
    scalarField dfdx1(nEqns());
    scalarSquareMatrix d2fdxdyAnalytic(nEqns());
    jacobian(x, y, li, dfdx1, d2fdxdyAnalytic);

    // Compare derivatives
    Info<< "[derivatives] dfdx = " << dfdx0 << nl;
    Info<< "[   jacobian] dfdx = " << dfdx1 << nl;

    // Construct a Jacobian using the finite differences and the derivatives
    // method
    scalarSquareMatrix d2fdxdyFiniteDifference(nEqns());
    scalarField y0(y), y1(y);
    for (label i = 0; i < nEqns(); ++ i)
    {
        y0[i] = y[i] - dy[i];
        y1[i] = y[i] + dy[i];

        derivatives(x, y0, li, dfdx0);
        derivatives(x, y1, li, dfdx1);

        for (label j = 0; j < nEqns(); j++)
        {
            d2fdxdyFiniteDifference(j, i) = (dfdx1[j] - dfdx0[j])/(2*dy[i]);
        }

        y0[i] = y[i];
        y1[i] = y[i];
    }

    for (label i = 0; i < nEqns(); ++ i)
    {
        Info<< "[derivatives] d2fdxdy[" << i << "] = "
            << UList<scalar>(d2fdxdyFiniteDifference[i], nEqns()) << nl;
        Info<< "[   jacobian] d2fdxdy[" << i << "] = "
            << UList<scalar>(d2fdxdyAnalytic[i], nEqns()) << nl;
    }
}


// ************************************************************************* //
