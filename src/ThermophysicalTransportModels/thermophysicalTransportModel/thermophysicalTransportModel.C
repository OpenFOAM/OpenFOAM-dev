/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "thermophysicalTransportModel.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermophysicalTransportModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermophysicalTransportModel::thermophysicalTransportModel
(
    const compressibleMomentumTransportModel& momentumTransport
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName
            (
                typeName, momentumTransport.alphaRhoPhi().group()
            ),
            momentumTransport.time().constant(),
            momentumTransport.mesh(),
            //***HGW IOobject::MUST_READ_IF_MODIFIED,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),

    momentumTransportModel_(momentumTransport)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::thermophysicalTransportModel::read()
{
    return regIOobject::read();
}


void Foam::thermophysicalTransportModel::correct()
{}


// ************************************************************************* //
