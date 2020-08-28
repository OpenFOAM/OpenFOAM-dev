/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::thermocapillaryForce

Description
    Thermocapillary force

SourceFiles
    thermocapillaryForce.C

\*---------------------------------------------------------------------------*/

#ifndef thermocapillaryForce_H
#define thermocapillaryForce_H

#include "force.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class thermocapillaryForce Declaration
\*---------------------------------------------------------------------------*/

class thermocapillaryForce
:
    public force
{
public:

    //- Runtime type information
    TypeName("thermocapillary");


    // Constructors

        //- Construct from surface film model
        thermocapillaryForce
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        thermocapillaryForce(const thermocapillaryForce&) = delete;


    //- Destructor
    virtual ~thermocapillaryForce();


    // Member Functions

        // Evolution

            //- Correct
            virtual tmp<fvVectorMatrix> correct(volVectorField& U);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const thermocapillaryForce&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
