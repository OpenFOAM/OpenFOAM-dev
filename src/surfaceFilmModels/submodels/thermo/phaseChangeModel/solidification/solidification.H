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
    Foam::surfaceFilmModels::solidification

Description
    Solidification phase change model where all film mass is converted when the
    local temperature > activation temperature.  The latent heat is
    assumed to be removed by heat-transfer to the wall.

SourceFiles
    solidification.C

\*---------------------------------------------------------------------------*/

#ifndef solidification_H
#define solidification_H

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                       Class solidification Declaration
\*---------------------------------------------------------------------------*/

class solidification
:
    public phaseChangeModel
{
protected:

    // Protected data

        //- Temperature at which solidification starts
        scalar T0_;

        //- Solidification limiter
        //  Maximum fraction of film which can solidify in a time-step
        scalar maxSolidificationFrac_;

        //- Solidification limiter
        //  Maximum rate at which the film can solidify
        dimensionedScalar maxSolidificationRate_;

        //- Accumulated solid mass [kg]
        volScalarField::Internal mass_;

        //- Accumulated solid thickness [m]
        volScalarField::Internal thickness_;


public:

    //- Runtime type information
    TypeName("solidification");


    // Constructors

        //- Construct from surface film model
        solidification(surfaceFilm& film, const dictionary& dict);

        //- Disallow default bitwise copy construction
        solidification(const solidification&) = delete;


    //- Destructor
    virtual ~solidification();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correctModel
            (
                const scalar dt,
                scalarField& availableMass,
                scalarField& dMass,
                scalarField& dEnergy
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solidification&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
