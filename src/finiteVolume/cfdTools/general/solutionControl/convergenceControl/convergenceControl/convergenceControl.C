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

#include "convergenceControl.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(convergenceControl, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::convergenceControl::getInitialResiduals
(
    const fvMesh& mesh,
    const word& fieldName,
    const label solvei,
    ITstream& data,
    scalar& r0,
    scalar& r
)
{
    getInitialTypeResiduals<scalar>(mesh, fieldName, solvei, data, r0, r);
    getInitialTypeResiduals<vector>(mesh, fieldName, solvei, data, r0, r);
    getInitialTypeResiduals<sphericalTensor>
    (
        mesh,
        fieldName,
        solvei,
        data,
        r0,
        r
    );
    getInitialTypeResiduals<symmTensor>(mesh, fieldName, solvei, data, r0, r);
    getInitialTypeResiduals<tensor>(mesh, fieldName, solvei, data, r0, r);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::convergenceControl::convergenceControl(const solutionControl& control)
:
    control_(control)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::convergenceControl::~convergenceControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::convergenceControl::converged()
{
    if
    (
        control_.time().timeIndex() != control_.time().startTimeIndex()
     && criteriaSatisfied()
    )
    {
        Info<< nl << control_.algorithmName() << " solution converged in "
            << control_.time().timeName() << " iterations" << nl << endl;

        return true;
    }

    return false;
}


bool Foam::convergenceControl::endIfConverged(Time& time)
{
    if (converged())
    {
        if (time.writeTime())
        {
            time.stopAt(Time::saNoWriteNow);
            time.setEndTime(time);
        }
        else
        {
            time.writeAndEnd();
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
