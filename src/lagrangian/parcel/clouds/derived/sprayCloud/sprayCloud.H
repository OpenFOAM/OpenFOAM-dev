/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Class
    Foam::sprayCloud

Description
    Cloud class to introduce reacting spray parcels

\*---------------------------------------------------------------------------*/

#ifndef sprayCloud_H
#define sprayCloud_H

#include "ParcelCloudBase.H"
#include "MomentumCloud.H"
#include "ThermoCloud.H"
#include "ReactingCloud.H"
#include "SprayCloud.H"
#include "ParcelCloud.H"
#include "sprayParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef ParcelCloud
    <
        SprayCloud
        <
            ReactingCloud
            <
                ThermoCloud
                <
                    MomentumCloud
                    <
                        ParcelCloudBase
                        <
                            sprayParcel
                        >
                    >
                >
            >
        >
    > sprayCloud;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
