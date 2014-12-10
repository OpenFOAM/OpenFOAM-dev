/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "pitchForkRing.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace tetherPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pitchForkRing, 0);

addToRunTimeSelectionTable
(
    tetherPotential,
    pitchForkRing,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pitchForkRing::pitchForkRing
(
    const word& name,
    const dictionary& tetherPotentialProperties
)
:
    tetherPotential(name, tetherPotentialProperties),
    pitchForkRingCoeffs_
    (
        tetherPotentialProperties.subDict(typeName + "Coeffs")
    ),
    mu_(readScalar(pitchForkRingCoeffs_.lookup("mu"))),
    alpha_(readScalar(pitchForkRingCoeffs_.lookup("alpha"))),
    rOrbit_(readScalar(pitchForkRingCoeffs_.lookup("rOrbit")))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar pitchForkRing::energy(const vector r) const
{
    scalar p = sqrt(r.x()*r.x() + r.y()*r.y());

    scalar pMinusRSqr = sqr(p - rOrbit_);

    return
       -0.5 * mu_ * pMinusRSqr
      + 0.25 * pMinusRSqr * pMinusRSqr
      + 0.5 * alpha_ * r.z() * r.z();
}


vector pitchForkRing::force(const vector r) const
{
    scalar p = sqrt(r.x()*r.x() + r.y()*r.y());

    scalar pMinusR = (p - rOrbit_);

    return vector
    (
        (mu_ - sqr(pMinusR)) * pMinusR * r.x()/(p + VSMALL),
        (mu_ - sqr(pMinusR)) * pMinusR * r.y()/(p + VSMALL),
      - alpha_ * r.z()
    );
}


bool pitchForkRing::read(const dictionary& tetherPotentialProperties)
{
    tetherPotential::read(tetherPotentialProperties);

    pitchForkRingCoeffs_ =
        tetherPotentialProperties.subDict(typeName + "Coeffs");

    pitchForkRingCoeffs_.lookup("mu") >> mu_;
    pitchForkRingCoeffs_.lookup("alpha") >> alpha_;
    pitchForkRingCoeffs_.lookup("rOrbit") >> rOrbit_;

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace tetherPotentials
} // End namespace Foam

// ************************************************************************* //
