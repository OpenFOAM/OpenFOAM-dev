/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "searchableSurfaceFeatures.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(searchableSurfaceFeatures, 0);
defineRunTimeSelectionTable(searchableSurfaceFeatures, dict);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::searchableSurfaceFeatures>
Foam::searchableSurfaceFeatures::New
(
    const searchableSurface& surface,
    const dictionary& dict
)
{
    word searchableSurfaceFeaturesType = surface.type() + "Features";

    dictConstructorTable::iterator cstrIter =
        dictConstructorTablePtr_->find(searchableSurfaceFeaturesType);

    if (cstrIter == dictConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown searchableSurfaceFeatures type "
            << searchableSurfaceFeaturesType << endl << endl
            << "Valid searchableSurfaceFeatures types : " << endl
            << dictConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<searchableSurfaceFeatures>(cstrIter()(surface, dict));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaceFeatures::searchableSurfaceFeatures
(
    const searchableSurface& surface,
    const dictionary& dict
)
:
    surface_(surface),
    dict_(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableSurfaceFeatures::~searchableSurfaceFeatures()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
