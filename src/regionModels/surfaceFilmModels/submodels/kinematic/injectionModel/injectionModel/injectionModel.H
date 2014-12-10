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
    Foam::regionModels::surfaceFilmModels::injectionModel

Description
    Base class for film injection models, handling mass transfer from the
    film.

SourceFiles
    injectionModel.C
    injectionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef injectionModel_H
#define injectionModel_H

#include "filmSubModelBase.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class injectionModel Declaration
\*---------------------------------------------------------------------------*/

class injectionModel
:
    public filmSubModelBase
{
private:

    // Private data

        //- Injected mass
        scalar injectedMass_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        injectionModel(const injectionModel&);

        //- Disallow default bitwise assignment
        void operator=(const injectionModel&);


protected:

    // Protected Member Functions

        //- Add to injected mass
        void addToInjectedMass(const scalar dMass);

        //- Correct
        void correct();


public:

    //- Runtime type information
    TypeName("injectionModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             injectionModel,
             dictionary,
             (
                surfaceFilmModel& owner,
                const dictionary& dict
             ),
             (owner, dict)
         );

    // Constructors

        //- Construct null
        injectionModel(surfaceFilmModel& owner);

        //- Construct from type name, dictionary and surface film model
        injectionModel
        (
            const word& modelType,
            surfaceFilmModel& owner,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected injection model
        static autoPtr<injectionModel> New
        (
            surfaceFilmModel& owner,
            const dictionary& dict,
            const word& mdoelType
        );


    //- Destructor
    virtual ~injectionModel();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct
            (
                scalarField& availableMass,
                scalarField& massToInject,
                scalarField& diameterToInject
            ) = 0;

            //- Return the total mass injected
            scalar injectedMassTotal() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
