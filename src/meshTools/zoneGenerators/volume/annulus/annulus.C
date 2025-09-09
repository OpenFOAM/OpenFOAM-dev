/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "annulus.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(annulus, 0);
        addToRunTimeSelectionTable(zoneGenerator, annulus, dictionary);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline bool Foam::zoneGenerators::annulus::contains(const point& p) const
{
    const vector d = p - point1_;
    const scalar magda = d & axis_;

    if ((magda > 0) && (magda < magAxis2_))
    {
        const scalar d2 = magSqr(d) - sqr(magda)/magAxis2_;

        if ((d2 < outerRadius2_) && (d2 > innerRadius2_))
        {
            return true;
        }
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::annulus::annulus
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    volume(name, mesh, dict),
    point1_(dict.lookup<point>("point1", dimLength)),
    point2_(dict.lookup<point>("point2", dimLength)),
    outerRadius_(dict.lookup<scalar>("outerRadius", dimLength)),
    innerRadius_(dict.lookup<scalar>("innerRadius", dimLength)),
    axis_(point2_ - point1_),
    magAxis2_(magSqr(axis_)),
    outerRadius2_(sqr(outerRadius_)),
    innerRadius2_(sqr(innerRadius_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::annulus::~annulus()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::annulus::generate() const
{
    return volume::generate(*this);
}


// ************************************************************************* //
