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
    Foam::regionModels::surfaceFilmModels::filmMomentumTransportModel

Description
    Base class for film turbulence models

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
                 Class filmMomentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

class filmMomentumTransportModel
:
    public filmSubModelBase
{
public:

    //- Runtime type information
    TypeName("filmMomentumTransportModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             filmMomentumTransportModel,
             dictionary,
             (
                surfaceFilmRegionModel& film,
                const dictionary& dict
             ),
             (film, dict)
         );

    // Constructors

        //- Construct null
        filmMomentumTransportModel(surfaceFilmRegionModel& film);

        //- Construct from type name, dictionary and surface film model
        filmMomentumTransportModel
        (
            const word& modelType,
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        filmMomentumTransportModel(const filmMomentumTransportModel&) = delete;


    // Selectors

        //- Return a reference to the selected injection model
        static autoPtr<filmMomentumTransportModel> New
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~filmMomentumTransportModel();


    // Member Functions

        // Evolution

            //- Return the film surface velocity
            virtual tmp<volVectorField::Internal> Us() const = 0;

            //- Return the film turbulence viscosity
            virtual tmp<volScalarField> mut() const = 0;

            //- Correct/update the model
            virtual void correct() = 0;

            //- Return the source for the film momentum equation
            virtual tmp<fvVectorMatrix> Su(volVectorField& U) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const filmMomentumTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
