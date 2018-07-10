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
    Foam::regionModels::surfaceFilmModels::
        perturbedTemperatureDependentContactAngleForce

Description
    Temperature dependent contact angle force with a stochastic perturbation

    The contact angle in degrees is specified as a Foam::Function1 type, to
    enable the use of, e.g.  contant, polynomial, table values and the
    stochastic perturbation obtained from a
    Foam::distributionModels::distributionModel

See also
    Foam::regionModels::surfaceFilmModels::contactAngleForce
    Foam::regionModels::surfaceFilmModels::temperatureDependentContactAngleForce
    Foam::regionModels::surfaceFilmModels::distributionContactAngleForce
    Foam::Function1Types
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
namespace surfaceFilmModels
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


    // Private member functions

        //- Disallow default bitwise copy construct
        perturbedTemperatureDependentContactAngleForce
        (
            const perturbedTemperatureDependentContactAngleForce&
        );

        //- Disallow default bitwise assignment
        void operator=(const perturbedTemperatureDependentContactAngleForce&);


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
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~perturbedTemperatureDependentContactAngleForce();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
