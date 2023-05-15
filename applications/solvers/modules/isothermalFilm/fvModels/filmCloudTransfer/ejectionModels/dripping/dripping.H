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
    Foam::filmEjectionModels::dripping

Description
    Dripping film to cloud ejection transfer model

    On an inverted surface if the film thickness is sufficient to generate a
    valid parcel the equivalent mass is removed from the film and transfered to
    the cloud as a parcel containing droplets with a diameter obtained from
    the specified parcelDistribution.

Usage
    Example usage:
    \verbatim
    filmCloudTransfer
    {
        type    filmCloudTransfer;

        libs    ("libfilmCloudTransfer.so");

        ejection
        {
            model   dripping;

            deltaStable 5e-4;

            minParticlesPerParcel 10;

            parcelDistribution
            {
                type            RosinRammler;
                Q               0;
                min             1e-3;
                max             2e-3;
                d               7.5e-05;
                n               0.5;
            }
        }
    }
    \endverbatim

SourceFiles
    dripping.C

\*---------------------------------------------------------------------------*/

#ifndef dripping_H
#define dripping_H

#include "ejectionModel.H"
#include "distribution.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace filmEjectionModels
{

/*---------------------------------------------------------------------------*\
                      Class dripping Declaration
\*---------------------------------------------------------------------------*/

class dripping
:
    public ejectionModel
{
    // Private Member Data

        //- Stable film thickness - drips only formed if thickness
        //  exceeds this threshold value
        scalar deltaStable_;

        //- Minimum number of droplets per parcel
        scalar minParticlesPerParcel_;

        //- Random number generator
        Random rndGen_;

        //- Parcel size PDF model
        const autoPtr<distribution> parcelDistribution_;


public:

    //- Runtime type information
    TypeName("dripping");


    // Constructors

        //- Construct from dictionary and film model
        dripping
        (
            const dictionary& dict,
            const solvers::isothermalFilm& film
        );

        //- Disallow default bitwise copy construction
        dripping(const dripping&) = delete;


    //- Destructor
    virtual ~dripping();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dripping&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace filmEjectionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
