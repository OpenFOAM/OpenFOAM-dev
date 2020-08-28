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
    Foam::regionModels::surfaceFilmModels::forceList

Description
    List container for film sources

SourceFiles
    forceList.C

\*---------------------------------------------------------------------------*/

#ifndef forceList_H
#define forceList_H

#include "PtrList.H"
#include "force.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                          Class forceList Declaration
\*---------------------------------------------------------------------------*/

class forceList
:
    public PtrList<force>
{
public:

    // Constructors

        //- Construct null
        forceList(surfaceFilmRegionModel& film);

        //- Construct from type name, dictionary and surface film model
        forceList
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~forceList();


    // Member Functions

        //- Return (net) force system
        tmp<fvVectorMatrix> correct(volVectorField& U);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
