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
    Foam::solidification

Description
    Solidification phase change model where all film mass is converted when
    the local temperature > activation temperature.

SourceFiles
    solidification.C

\*---------------------------------------------------------------------------*/

#ifndef solidification_H
#define solidification_H

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
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
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        solidification(const solidification&);

        //- Disallow default bitwise assignment
        void operator=(const solidification&);


protected:

    // Protected data

        //- Temperature at which solidification starts
        scalar T0_;

        //- Latent heat of solidification [J/kg]
        scalar L_;

        //- Under-relaxation parameter for solidification process (0-1)
        scalar alpha_;

        //- Accumulated solid mass [kg]
        volScalarField mass_;

        //- Accumulated solid thickness [m]
        volScalarField thickness_;


public:

    //- Runtime type information
    TypeName("solidification");


    // Constructors

        //- Construct from surface film model
        solidification(surfaceFilmModel& owner, const dictionary& dict);


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
