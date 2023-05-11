/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::surfaceFilmModels::drippingEjection

Description
    Film Dripping mass transfer model.

    If the film mass exceeds that needed to generate a valid parcel, the
    equivalent mass is removed from the film.

    New parcel diameters are sampled from a PDF.

SourceFiles
    drippingEjection.C

\*---------------------------------------------------------------------------*/

#ifndef drippingEjection_H
#define drippingEjection_H

#include "ejectionModel.H"
#include "distribution.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class drippingEjection Declaration
\*---------------------------------------------------------------------------*/

class drippingEjection
:
    public ejectionModel
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
        const autoPtr<distribution> parcelDistribution_;

        //- Diameters of particles to eject into the dripping
        scalarList diameter_;


public:

    //- Runtime type information
    TypeName("drippingEjection");


    // Constructors

        //- Construct from surface film model
        drippingEjection
        (
            surfaceFilm& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        drippingEjection(const drippingEjection&) = delete;


    //- Destructor
    virtual ~drippingEjection();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct
            (
                scalarField& availableMass,
                scalarField& massToEject,
                scalarField& diameterToEject
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const drippingEjection&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
