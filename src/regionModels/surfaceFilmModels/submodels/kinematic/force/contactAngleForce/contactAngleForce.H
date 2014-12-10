/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::contactAngleForce

Description
    Film contact angle force

    The effect of the contact angle force can be ignored over a specified
    distance from patches.

SourceFiles
    contactAngleForce.C

\*---------------------------------------------------------------------------*/

#ifndef contactAngleForce_H
#define contactAngleForce_H

#include "force.H"
#include "distributionModel.H"
#include "cachedRandom.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class contactAngleForce Declaration
\*---------------------------------------------------------------------------*/

class contactAngleForce
:
    public force
{
private:

    // Private Data

        //- Coefficient applied to the contact angle force
        scalar Ccf_;

        //- Random number generator
        cachedRandom rndGen_;

        //- Parcel size PDF model
        const autoPtr<distributionModels::distributionModel> distribution_;

        //- Mask over which force is applied
        volScalarField mask_;


    // Private member functions

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        contactAngleForce(const contactAngleForce&);

        //- Disallow default bitwise assignment
        void operator=(const contactAngleForce&);


public:

    //- Runtime type information
    TypeName("contactAngle");


    // Constructors

        //- Construct from surface film model
        contactAngleForce
        (
            surfaceFilmModel& owner,
            const dictionary& dict
        );


    //- Destructor
    virtual ~contactAngleForce();


    // Member Functions

        // Evolution

            //- Correct
            virtual tmp<fvVectorMatrix> correct(volVectorField& U);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
