/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::ejectionModelList

Description
    List container for film ejection models

SourceFiles
    ejectionModelList.C

\*---------------------------------------------------------------------------*/

#ifndef ejectionModelList_H
#define ejectionModelList_H

#include "PtrList.H"
#include "ejectionModel.H"
#include "filmSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                     Class ejectionModelList Declaration
\*---------------------------------------------------------------------------*/

class ejectionModelList
:
    public PtrList<ejectionModel>,
    public filmSubModelBase
{
    // Private Data

        //- List of mass ejected per patch
        scalarField massEjected_;


public:

    // Constructors

        //- Construct null
        ejectionModelList(surfaceFilmRegionModel& film);

        //- Construct from type name, dictionary and surface film model
        ejectionModelList
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        ejectionModelList(const ejectionModelList&) = delete;


    //- Destructor
    virtual ~ejectionModelList();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct
            (
                scalarField& availableMass,
                volScalarField& massToEject,
                volScalarField& diameterToEject
            );


        // I-O

            //- Provide some info
            virtual void info(Ostream& os);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ejectionModelList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
