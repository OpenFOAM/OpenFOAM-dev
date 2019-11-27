/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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

#include "piecewiseLinearRamp.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(piecewiseLinearRamp, 0);
addToRunTimeSelectionTable
(
    faceAreaWeightModel,
    piecewiseLinearRamp,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

piecewiseLinearRamp::piecewiseLinearRamp
(
    const dictionary& faceAreaWeightDict
)
:
    faceAreaWeightModel(typeName, faceAreaWeightDict),
    lAF_(coeffDict().lookup<scalar>("lowerAreaFraction")),
    uAF_(coeffDict().lookup<scalar>("upperAreaFraction"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar piecewiseLinearRamp::faceAreaWeight(scalar faceAreaFraction) const
{
    if (faceAreaFraction < lAF_)
    {
        return 0;
    }
    else if (faceAreaFraction < uAF_)
    {
        return faceAreaFraction/((uAF_ - lAF_)) - lAF_/(uAF_ - lAF_);
    }
    else
    {
        return 1;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
