/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::BrunDrippingInjection

Description
    Film Dripping mass transfer model.

    If the film thickness exceeds the critical value needed to generate a valid
    parcel, the equivalent mass is removed from the film.  The critical film
    thickness is calculated from the Rayleigh-Taylor stability analysis
    of film flow on an inclined plane by Brun et.al.

    Reference:
    \verbatim
        Brun, P. T., Damiano, A., Rieu, P., Balestra, G., & Gallaire, F. (2015).
        Rayleigh-Taylor instability under an inclined plane.
        Physics of Fluids (1994-present), 27(8), 084107.
    \endverbatim

    New parcel diameters are sampled from a PDF.

SourceFiles
    BrunDrippingInjection.C

\*---------------------------------------------------------------------------*/

#ifndef BrunDrippingInjection_H
#define BrunDrippingInjection_H

#include "injectionModel.H"
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
                    Class BrunDrippingInjection Declaration
\*---------------------------------------------------------------------------*/

class BrunDrippingInjection
:
    public injectionModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        BrunDrippingInjection(const BrunDrippingInjection&);

        //- Disallow default bitwise assignment
        void operator=(const BrunDrippingInjection&);


protected:

    // Protected data

        //- Critical non-dimensional interface velocity
        //  Coefficient in the film angle stability function
        // defaults to 1.62208
        scalar ubarStar_;

        //- Stable film thickness - drips only formed if thickness
        //  execeeds this threhold value
        scalar deltaStable_;

        //- Number of particles per parcel
        scalar particlesPerParcel_;

        //- Random number generator
        cachedRandom rndGen_;

        //- Parcel size PDF model
        const autoPtr<distributionModels::distributionModel>
            parcelDistribution_;

        //- Diameters of particles to inject into the dripping
        scalarList diameter_;


public:

    //- Runtime type information
    TypeName("BrunDrippingInjection");


    // Constructors

        //- Construct from surface film model
        BrunDrippingInjection
        (
            surfaceFilmModel& owner,
            const dictionary& dict
        );


    //- Destructor
    virtual ~BrunDrippingInjection();


    // Member Functions

        //- Correct
        virtual void correct
        (
            scalarField& availableMass,
            scalarField& massToInject,
            scalarField& diameterToInject
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
