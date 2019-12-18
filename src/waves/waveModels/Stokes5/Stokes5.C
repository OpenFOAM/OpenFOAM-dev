/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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

#include "Stokes5.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(Stokes5, 0);
    addToRunTimeSelectionTable(waveModel, Stokes5, objectRegistry);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::Stokes5::Stokes5
(
    const objectRegistry& db,
    const dictionary& dict
)
:
    Stokes2(db, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveModels::Stokes5::~Stokes5()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::waveModels::Stokes5::elevation
(
    const scalar t,
    const scalarField& x
) const
{
    static const scalar kdGreat = log(great);
    const scalar kd = min(max(k()*depth(), - kdGreat), kdGreat);
    const scalar ka = k()*amplitude(t);

    const scalar S = deep() ? 0 : 1/cosh(2*kd), T = deep() ? 1 : tanh(kd);

    const scalar B31 =
      - 3.0/8/pow3(1 - S)
       *(
            1 + 3*S + 3*sqr(S) + 2*pow3(S)
        );

    const scalar B42 =
        1.0/6/T/(3 + 2*S)/pow4(1 - S)
       *(
            6 - 26*S - 182*sqr(S) - 204*pow3(S) - 25*pow4(S) + 26*pow5(S)
        );

    const scalar B44 =
        1.0/24/T/(3 + 2*S)/pow4(1 - S)
       *(
            24 + 92*S + 122*sqr(S) + 66*pow3(S) + 67*pow4(S) + 34*pow5(S)
        );

    const scalar B53 =
        9.0/128/(3 + 2*S)/(4 + S)/pow6(1 - S)
       *(
            132 + 17*S - 2216*sqr(S) - 5897*pow3(S) - 6292*pow4(S)
          - 2687*pow5(S) + 194*pow6(S) + 467*S*pow6(S) + 82*sqr(pow4(S))
        );

    const scalar B55 =
        5.0/384/(3 + 2*S)/(4 + S)/pow6(1 - S)
       *(
            300 + 1579*S + 3176*sqr(S) + 2949*pow3(S) + 1188*pow4(S)
          + 675*pow5(S) + 1326*pow6(S) + 827*S*pow6(S) + 130*sqr(pow4(S))
        );

    if (debug)
    {
        Info<< "B31 = " << B31 << endl
            << "B42 = " << B42 << endl
            << "B44 = " << B44 << endl
            << "B53 = " << B53 << endl
            << "B55 = " << B55 << endl;
    }

    const scalarField phi(angle(t, x));

    return
        Stokes2::elevation(t, x)
      + (1/k())
       *(
           pow3(ka)*B31*(cos(phi) - cos(3*phi))
         + pow4(ka)*(B42*cos(2*phi) + B44*cos(4*phi))
         + pow5(ka)*(- (B53 + B55)*cos(phi) + B53*cos(3*phi) + B55*cos(5*phi))
        );
}


Foam::tmp<Foam::vector2DField> Foam::waveModels::Stokes5::velocity
(
    const scalar t,
    const vector2DField& xz
) const
{
    static const scalar kdGreat = log(great);
    const scalar kd = min(max(k()*depth(), - kdGreat), kdGreat);
    const scalar ka = k()*amplitude(t);

    const scalar S = deep() ? 0 : 1/cosh(2*kd);
    const scalar SByA11 = deep() ? 0 : S*sinh(kd);

    const scalar A31ByA11 =
        1.0/8/pow3(1 - S)
       *(
          - 4 - 20*S + 10*sqr(S) - 13*pow3(S)
        );

    const scalar A33ByA11 =
        1.0/8/pow3(1 - S)
       *(
          - 2*sqr(S) + 11*pow3(S)
        );

    const scalar A42ByA11 =
        SByA11/24/pow5(1 - S)
       *(
            12 - 14*S - 264*sqr(S) - 45*pow3(S) - 13*pow4(S)
        );

    const scalar A44ByA11 =
        SByA11/48/(3 + 2*S)/pow5(1 - S)
       *(
            10*sqr(S) - 174*pow3(S) + 291*pow4(S) + 278*pow5(S)
        );

    const scalar A51ByA11 =
        1.0/64/(3 + 2*S)/(4 + S)/pow6(1 - S)
       *(
         - 1184 + 32*S + 13232*sqr(S) + 21712*pow3(S) + 20940*pow4(S)
         + 12554*pow5(S) - 500*pow6(S) - 3341*S*pow6(S) - 670*sqr(pow4(S))
        );

    const scalar A53ByA11 =
        1.0/32/(3 + 2*S)/pow6(1 - S)
       *(
            4*S + 105*sqr(S) + 198*pow3(S) - 1376*pow4(S) - 1302*pow5(S)
          - 117*pow6(S) + 58*S*pow6(S)
        );

    const scalar A55ByA11 =
        1.0/64/(3 + 2*S)/(4 + S)/pow6(1 - S)
       *(
          - 6*pow3(S) + 272*pow4(S) - 1552*pow5(S) + 852*pow6(S)
          + 2029*S*pow6(S) + 430*sqr(pow4(S))
        );

    if (debug)
    {
        const scalar A11 = 1/sinh(kd);
        Info<< "A31 = " << A31ByA11*A11 << endl
            << "A33 = " << A33ByA11*A11 << endl
            << "A42 = " << A42ByA11*A11 << endl
            << "A44 = " << A44ByA11*A11 << endl
            << "A51 = " << A51ByA11*A11 << endl
            << "A53 = " << A53ByA11*A11 << endl
            << "A55 = " << A55ByA11*A11 << endl;
    }

    const vector2DField v1(vi(1, t, xz)), v3(vi(3, t, xz));

    return
        Stokes2::velocity(t, xz)
      + celerity()
       *(
            pow3(ka)*(A31ByA11*v1 + A33ByA11*v3)
          + pow4(ka)*(A42ByA11*vi(2, t, xz) + A44ByA11*vi(4, t, xz))
          + pow5(ka)*(A51ByA11*v1 + A53ByA11*v3 + A55ByA11*vi(5, t, xz))
        );
}


// ************************************************************************* //
