/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::temperatureDependentContactAngleForce

Description
    Temperature dependent contact angle force

    The contact angle in degrees is specified as a Foam::Function1 type, to
    enable the use of, e.g.  contant, polynomial, table values.

See also
    Foam::regionModels::surfaceFilmModels::contactAngleForce
    Foam::Function1Types

SourceFiles
    temperatureDependentContactAngleForce.C

\*---------------------------------------------------------------------------*/

#ifndef temperatureDependentContactAngleForce_H
#define temperatureDependentContactAngleForce_H

#include "contactAngleForce.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
            Class temperatureDependentContactAngleForce Declaration
\*---------------------------------------------------------------------------*/

class temperatureDependentContactAngleForce
:
    public contactAngleForce
{
    // Private Data

        //- Contact angle function
        autoPtr<Function1<scalar>> thetaPtr_;


    // Private member functions

        //- Disallow default bitwise copy construct
        temperatureDependentContactAngleForce
        (
            const temperatureDependentContactAngleForce&
        );

        //- Disallow default bitwise assignment
        void operator=(const temperatureDependentContactAngleForce&);


protected:

        //- Return the contact angle field
        virtual tmp<volScalarField> theta() const;


public:

    //- Runtime type information
    TypeName("temperatureDependentContactAngle");


    // Constructors

        //- Construct from surface film model
        temperatureDependentContactAngleForce
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~temperatureDependentContactAngleForce();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
