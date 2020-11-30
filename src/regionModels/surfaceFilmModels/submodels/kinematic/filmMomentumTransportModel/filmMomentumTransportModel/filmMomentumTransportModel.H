/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::momentumTransportModel

Description
    Base class for film momentum transport models

SourceFiles
    filmMomentumTransportModel.C
    filmMomentumTransportModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef filmMomentumTransportModel_H
#define filmMomentumTransportModel_H

#include "filmSubModelBase.H"
#include "runTimeSelectionTables.H"
#include "fvMatricesFwd.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                 Class momentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

class momentumTransportModel
:
    public filmSubModelBase
{
public:

    //- Runtime type information
    TypeName("momentumTransport");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             momentumTransportModel,
             dictionary,
             (
                surfaceFilmRegionModel& film,
                const dictionary& dict
             ),
             (film, dict)
         );

    // Constructors

        //- Construct null
        momentumTransportModel(surfaceFilmRegionModel& film);

        //- Construct from type name, dictionary and surface film model
        momentumTransportModel
        (
            const word& modelType,
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        momentumTransportModel(const momentumTransportModel&) = delete;


    // Selectors

        //- Return a reference to the selected injection model
        static autoPtr<momentumTransportModel> New
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~momentumTransportModel();


    // Member Functions

        // Evolution

            //- Return the film surface velocity
            virtual tmp<volVectorField::Internal> Us() const = 0;

            //- Correct/update the model
            virtual void correct() = 0;

            //- Return the source for the film momentum equation
            virtual tmp<fvVectorMatrix> Su(volVectorField& U) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const momentumTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
