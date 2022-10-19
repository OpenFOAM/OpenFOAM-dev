/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "filmSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmSubModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

filmSubModelBase::filmSubModelBase(surfaceFilm& film)
:
    subModelBase(film.outputProperties()),
    filmModel_(film)
{}


filmSubModelBase::filmSubModelBase
(
    surfaceFilm& film,
    const dictionary& dict,
    const word& baseName,
    const word& modelType,
    const word& dictExt
)
:
    subModelBase
    (
        film.outputProperties(),
        dict,
        baseName,
        modelType,
        dictExt
    ),
    filmModel_(film)
{}


filmSubModelBase::filmSubModelBase
(
    const word& modelName,
    surfaceFilm& film,
    const dictionary& dict,
    const word& baseName,
    const word& modelType
)
:
    subModelBase
    (
        modelName,
        film.outputProperties(),
        dict,
        baseName,
        modelType
    ),
    filmModel_(film)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

filmSubModelBase::~filmSubModelBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool filmSubModelBase::writeTime() const
{
    return filmModel_.time().writeTime();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmSubModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
