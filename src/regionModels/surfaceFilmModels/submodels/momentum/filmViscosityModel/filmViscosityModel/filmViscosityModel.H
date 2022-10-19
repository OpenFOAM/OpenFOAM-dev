/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmSubModels::viscosityModel

Description
    Base class for surface film viscosity models

SourceFiles
    filmViscosityModel.C
    filmViscosityModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef filmViscosityModel_H
#define filmViscosityModel_H

#include "filmSubModelBase.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmSubModels
{

/*---------------------------------------------------------------------------*\
                     Class viscosityModel Declaration
\*---------------------------------------------------------------------------*/

class viscosityModel
:
    public filmSubModelBase
{
protected:

    // Protected Member Data

        //- Reference to the viscosity field
        volScalarField& mu_;


public:

    //- Runtime type information
    TypeName("viscosity");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             viscosityModel,
             dictionary,
             (
                surfaceFilm& film,
                const dictionary& dict,
                volScalarField& mu
             ),
             (film, dict, mu)
         );

    // Constructors

        //- Construct from type name, dictionary and surface film model
        viscosityModel
        (
            const word& modelType,
            surfaceFilm& film,
            const dictionary& dict,
            volScalarField& mu
        );

        //- Disallow default bitwise copy construction
        viscosityModel(const viscosityModel&) = delete;


    // Selectors

        //- Return a reference to the selected phase change model
        static autoPtr<viscosityModel> New
        (
            surfaceFilm& film,
            const dictionary& dict,
            volScalarField& mu
        );


    //- Destructor
    virtual ~viscosityModel();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct
            (
                const volScalarField& p,
                const volScalarField& T
            ) = 0;


        // I-O

            //- Provide some feedback
            virtual void info(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const viscosityModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmSubModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
