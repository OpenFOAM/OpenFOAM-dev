/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::force

Description
    Base class for film (stress-based) force models

SourceFiles
    force.C
    forceNew.C

\*---------------------------------------------------------------------------*/

#ifndef force_H
#define force_H

#include "filmSubModelBase.H"
#include "runTimeSelectionTables.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                          Class force Declaration
\*---------------------------------------------------------------------------*/

class force
:
    public filmSubModelBase
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        force(const force&);

        //- Disallow default bitwise assignment
        void operator=(const force&);


public:

    //- Runtime type information
    TypeName("force");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             force,
             dictionary,
             (
                surfaceFilmModel& owner,
                const dictionary& dict
             ),
             (owner, dict)
         );

    // Constructors

        //- Construct null
        force(surfaceFilmModel& owner);

        //- Construct from type name, dictionary and surface film model
        force
        (
            const word& modelType,
            surfaceFilmModel& owner,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected force model
        static autoPtr<force> New
        (
            surfaceFilmModel& owner,
            const dictionary& dict,
            const word& modelType
        );


    //- Destructor
    virtual ~force();


    // Member Functions

        // Evolution

            //- Correct
            virtual tmp<fvVectorMatrix> correct(volVectorField& U) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
