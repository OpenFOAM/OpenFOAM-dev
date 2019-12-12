/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::filmRadiationModel

Description
    Base class for film radiation models

SourceFiles
    filmRadiationModel.C
    filmRadiationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef filmRadiationModel_H
#define filmRadiationModel_H

#include "filmSubModelBase.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class filmRadiationModel Declaration
\*---------------------------------------------------------------------------*/

class filmRadiationModel
:
    public filmSubModelBase
{
public:

    //- Runtime type information
    TypeName("radiationModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             filmRadiationModel,
             dictionary,
             (
                 surfaceFilmRegionModel& film,
                 const dictionary& dict
             ),
             (film, dict)
         );

    // Constructors

        //- Construct null
        filmRadiationModel(surfaceFilmRegionModel& film);

        //- Construct from type name, dictionary and surface film model
        filmRadiationModel
        (
            const word& modelType,
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        filmRadiationModel(const filmRadiationModel&) = delete;


    // Selectors

        //- Return a reference to the selected phase change model
        static autoPtr<filmRadiationModel> New
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~filmRadiationModel();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct() = 0;

            //- Return the radiation sensible enthalpy source
            virtual tmp<volScalarField::Internal> Shs() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const filmRadiationModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
