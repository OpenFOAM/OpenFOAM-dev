/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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
    Foam::filmEjectionModels::BrunDripping

Description
    Brun dripping film to cloud ejection transfer model

    If the film thickness exceeds the critical value needed to generate one or
    more drops, the equivalent mass is removed from the film.  The critical film
    thickness is calculated from the Rayleigh-Taylor stability analysis of film
    flow on an inclined plane by Brun et.al.

    Reference:
    \verbatim
        Brun, P. T., Damiano, A., Rieu, P., Balestra, G., & Gallaire, F. (2015).
        Rayleigh-Taylor instability under an inclined plane.
        Physics of Fluids (1994-present), 27(8), 084107.
    \endverbatim

    The diameter of the drops formed are obtained from the local capillary
    length multiplied by the \c dCoeff coefficient which defaults to 3.3.

    Reference:
    \verbatim
        Lefebvre, A. (1988).
        Atomisation and sprays
        (Vol. 1040, No. 2756). CRC press.
    \endverbatim

Usage
    Example usage:
    \verbatim
    filmCloudTransfer
    {
        type    filmCloudTransfer;

        ejection
        {
            model   BrunDripping;

            deltaStable 5e-4;
        }
    }
    \endverbatim

SourceFiles
    BrunDripping.C

\*---------------------------------------------------------------------------*/

#ifndef BrunDripping_H
#define BrunDripping_H

#include "ejectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace filmEjectionModels
{

/*---------------------------------------------------------------------------*\
                    Class BrunDripping Declaration
\*---------------------------------------------------------------------------*/

class BrunDripping
:
    public ejectionModel
{
    // Private Member Data

        //- Critical non-dimensional interface velocity
        //  Coefficient in the film angle stability function.
        //  Defaults to 1.62208
        scalar ubarStar_;

        //- Coefficient relating the diameter of the drops formed to
        //  the capillary length.
        //  Defaults to 3.3
        scalar dCoeff_;

        //- Stable film thickness - drips only formed if thickness
        //  exceeds this threshold value
        scalar deltaStable_;

        //- Minimum number of droplets per parcel
        scalar minParticlesPerParcel_;


public:

    //- Runtime type information
    TypeName("BrunDripping");


    // Constructors

        //- Construct from dictionary and film model
        BrunDripping
        (
            const dictionary& dict,
            const solvers::isothermalFilm& film
        );

        //- Disallow default bitwise copy construction
        BrunDripping(const BrunDripping&) = delete;


    //- Destructor
    virtual ~BrunDripping();


    // Member Functions

        //- Correct
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const BrunDripping&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace filmEjectionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
