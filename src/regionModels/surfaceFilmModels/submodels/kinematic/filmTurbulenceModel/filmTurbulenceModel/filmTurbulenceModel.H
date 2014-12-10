/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::filmTurbulenceModel

Description
    Base class for film turbulence models

SourceFiles
    filmTurbulenceModel.C
    filmTurbulenceModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef filmTurbulenceModel_H
#define filmTurbulenceModel_H

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
                       Class filmTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class filmTurbulenceModel
:
    public filmSubModelBase
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        filmTurbulenceModel(const filmTurbulenceModel&);

        //- Disallow default bitwise assignment
        void operator=(const filmTurbulenceModel&);


public:

    //- Runtime type information
    TypeName("filmTurbulenceModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             filmTurbulenceModel,
             dictionary,
             (
                surfaceFilmModel& owner,
                const dictionary& dict
             ),
             (owner, dict)
         );

    // Constructors

        //- Construct null
        filmTurbulenceModel(surfaceFilmModel& owner);

        //- Construct from type name, dictionary and surface film model
        filmTurbulenceModel
        (
            const word& modelType,
            surfaceFilmModel& owner,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected injection model
        static autoPtr<filmTurbulenceModel> New
        (
            surfaceFilmModel& owner,
            const dictionary& dict
        );


    //- Destructor
    virtual ~filmTurbulenceModel();


    // Member Functions

        // Evolution

            //- Return the film surface velocity
            virtual tmp<volVectorField> Us() const = 0;

            //- Return the film turbulence viscosity
            virtual tmp<volScalarField> mut() const = 0;

            //- Correct/update the model
            virtual void correct() = 0;

            //- Return the source for the film momentum equation
            virtual tmp<fvVectorMatrix> Su(volVectorField& U) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
