/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2022 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmSubModels::
        perturbedTemperatureDependentContactAngleForce

Description
    Temperature dependent contact angle force with a stochastic perturbation

    The contact angle in degrees is specified as a Foam::Function1 type, to
    enable the use of, e.g.  constant, polynomial, table values and the
    stochastic perturbation obtained from a
    Foam::distributionModels::distributionModel

See also
    Foam::regionModels::surfaceFilmSubModels::contactAngleForce
    Foam::regionModels::surfaceFilmSubModels::\\
    temperatureDependentContactAngleForce
    Foam::regionModels::surfaceFilmSubModels::distributionContactAngleForce
    Foam::Function1s
    Foam::distributionModel

SourceFiles
    perturbedTemperatureDependentContactAngleForce.C

\*---------------------------------------------------------------------------*/

#ifndef perturbedTemperatureDependentContactAngleForce_H
#define perturbedTemperatureDependentContactAngleForce_H

#include "contactAngleForce.H"
#include "Function1.H"
#include "distributionModel.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmSubModels
{

/*---------------------------------------------------------------------------*\
       Class perturbedTemperatureDependentContactAngleForce Declaration
\*---------------------------------------------------------------------------*/

class perturbedTemperatureDependentContactAngleForce
:
    public contactAngleForce
{
    // Private Data

        //- Contact angle function
        autoPtr<Function1<scalar>> thetaPtr_;

        //- Random number generator
        Random rndGen_;

        //- Parcel size PDF model
        const autoPtr<distributionModel> distribution_;


protected:

        //- Return the contact angle field
        virtual tmp<volScalarField> theta() const;


public:

    //- Runtime type information
    TypeName("perturbedTemperatureDependentContactAngle");


    // Constructors

        //- Construct from surface film model
        perturbedTemperatureDependentContactAngleForce
        (
            surfaceFilm& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        perturbedTemperatureDependentContactAngleForce
        (
            const perturbedTemperatureDependentContactAngleForce&
        ) = delete;


    //- Destructor
    virtual ~perturbedTemperatureDependentContactAngleForce();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const perturbedTemperatureDependentContactAngleForce&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmSubModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
