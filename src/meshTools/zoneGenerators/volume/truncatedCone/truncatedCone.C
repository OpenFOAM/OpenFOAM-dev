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

#include "truncatedCone.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(truncatedCone, 0);
        addToRunTimeSelectionTable(zoneGenerator, truncatedCone, dictionary);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline bool Foam::zoneGenerators::truncatedCone::contains(const point& p) const
{
    const vector d = p - point1_;
    const scalar magda = d & axis_;

    if ((magda > 0) && (magda < magAxis2_))
    {
        const scalar d2 = magSqr(d) - sqr(magda)/magAxis2_;
        const scalar rAtp =
            radius1_ + (radius2_ - radius1_)*magda/magSqr(axis_);

        if (d2 <= sqr(rAtp))
        {
            return true;
        }
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::truncatedCone::truncatedCone
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    volume(name, mesh, dict),
    point1_(dict.lookup<point>("point1", dimLength)),
    point2_(dict.lookup<point>("point2", dimLength)),
    radius1_(dict.lookup<scalar>("radius1", dimLength)),
    radius2_(dict.lookup<scalar>("radius2", dimLength)),
    axis_(point2_ - point1_),
    magAxis2_(magSqr(axis_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::truncatedCone::~truncatedCone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::truncatedCone::generate() const
{
    return volume::generate(*this);
}


// ************************************************************************* //
