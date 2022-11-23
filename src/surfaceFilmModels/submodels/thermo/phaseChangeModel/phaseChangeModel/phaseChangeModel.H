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

Class
    Foam::regionModels::surfaceFilmSubModels::phaseChangeModel

Description
    Base class for surface film phase change models

SourceFiles
    phaseChangeModel.C
    phaseChangeModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeModel_H
#define phaseChangeModel_H

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
                      Class phaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeModel
:
    public filmSubModelBase
{
protected:

    // Protected Member Functions

        //- Latest mass transfer due to phase change
        scalar latestMassPC_;

        //- Total mass transfer due to phase change
        scalar totalMassPC_;


public:

    //- Runtime type information
    TypeName("phaseChange");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             phaseChangeModel,
             dictionary,
             (
                surfaceFilm& film,
                const dictionary& dict
             ),
             (film, dict)
         );

    // Constructors

        //- Construct null
        phaseChangeModel(surfaceFilm& film);

        //- Construct from type name, dictionary and surface film model
        phaseChangeModel
        (
            const word& modelType,
            surfaceFilm& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        phaseChangeModel(const phaseChangeModel&) = delete;


    // Selectors

        //- Return a reference to the selected phase change model
        static autoPtr<phaseChangeModel> New
        (
            surfaceFilm& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~phaseChangeModel();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct
            (
                const scalar dt,
                scalarField& availableMass,
                volScalarField& dMass,
                volScalarField& dEnergy
            );

            //- Correct
            virtual void correctModel
            (
                const scalar dt,
                scalarField& availableMass,
                scalarField& dMass,
                scalarField& dEnergy
            ) = 0;


        // I-O

            //- Provide some feedback
            virtual void info(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phaseChangeModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmSubModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
