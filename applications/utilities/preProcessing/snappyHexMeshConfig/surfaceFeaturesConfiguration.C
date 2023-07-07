/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "surfaceFeaturesConfiguration.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfaceFeaturesConfiguration::writeSurfaces()
{
    beginList(os_, "surfaces");

    forAll(surfaces_, i)
    {
        os_ << indent << surfaces_[i].file() << endl;
    }

    endList(os_);
}


void Foam::surfaceFeaturesConfiguration::writeSubsetFeatures()
{
    beginDict(os_, "subsetFeatures");
    os_ << indent << "nonManifoldEdges yes;" << endl;
    os_ << indent << "openEdges        yes;" << endl;
    endDict(os_);
}


void Foam::surfaceFeaturesConfiguration::writeTrimFeatures()
{
    beginDict(os_, "trimFeatures");
    os_ << indent << "minElem         0;" << endl;
    os_ << indent << "minLen          0;" << endl;
    endDict(os_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFeaturesConfiguration::surfaceFeaturesConfiguration
(
    const fileName& name,
    const fileName& dir,
    const Time& time,
    const meshingSurfaceList& surfaces
)
:
    caseFileConfiguration(name, dir, time),
    surfaces_(surfaces)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceFeaturesConfiguration::~surfaceFeaturesConfiguration()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::surfaceFeaturesConfiguration::write()
{
    dict_.writeHeader(os_, word("dictionary"));

    writeSurfaces();

    os_ << "includedAngle   150;" << nl << endl;

    writeSubsetFeatures();
    writeTrimFeatures();

    os_ << "writeObj        yes;";

    dict_.writeEndDivider(os_);
}


// ************************************************************************* //
