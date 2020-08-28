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
    Foam::regionModels::surfaceFilmModels::drippingInjection

Description
    Film Dripping mass transfer model.

    If the film mass exceeds that needed to generate a valid parcel, the
    equivalent mass is removed from the film.

    New parcel diameters are sampled from a PDF.

SourceFiles
    drippingInjection.C

\*---------------------------------------------------------------------------*/

#ifndef drippingInjection_H
#define drippingInjection_H

#include "injectionModel.H"
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
                      Class drippingInjection Declaration
\*---------------------------------------------------------------------------*/

class drippingInjection
:
    public injectionModel
{
protected:

    // Protected data

        //- Stable film thickness - drips only formed if thickness
        //  exceeds this threshold value
        scalar deltaStable_;

        //- Number of particles per parcel
        scalar particlesPerParcel_;

        //- Random number generator
        Random rndGen_;

        //- Parcel size PDF model
        const autoPtr<distributionModel>
            parcelDistribution_;

        //- Diameters of particles to inject into the dripping
        scalarList diameter_;


public:

    //- Runtime type information
    TypeName("drippingInjection");


    // Constructors

        //- Construct from surface film model
        drippingInjection
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        drippingInjection(const drippingInjection&) = delete;


    //- Destructor
    virtual ~drippingInjection();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct
            (
                scalarField& availableMass,
                scalarField& massToInject,
                scalarField& diameterToInject
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const drippingInjection&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
