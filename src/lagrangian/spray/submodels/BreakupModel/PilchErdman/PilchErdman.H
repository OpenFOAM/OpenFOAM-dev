/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::PilchErdman

Description
    Particle secondary breakup model, based on the reference:

    @verbatim
    Pilch, M. and Erdman, C.A.
    "Use of breakup time data and velocity history data
     to predict the maximum size of stable fragments for acceleration
     induced breakup of a liquid drop."
    Int. J. Multiphase Flows 13 (1987), 741-757
    @endverbatim

    The droplet fragment velocity is described by the equation:

    \f[
        V_d = V sqrt(epsilon)(B1 T + B2 T^2)
    \f]

    Where:
        V_d     : fragment velocity
        V       : magnitude of the relative velocity
        epsilon : density ratio (rho_carrier/rho_droplet)
        T       : characteristic break-up time
        B1, B2  : model input coefficients

    The authors suggest that:
        compressible flow   : B1 = 0.75*1.0; B2 = 3*0.116
        incompressible flow : B1 = 0.75*0.5; B2 = 3*0.0758


\*---------------------------------------------------------------------------*/

#ifndef PilchErdman_H
#define PilchErdman_H

#include "BreakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class PilchErdman Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PilchErdman
:
    public BreakupModel<CloudType>
{
private:

    // Private data

        scalar B1_;
        scalar B2_;


public:

    //- Runtime type information
    TypeName("PilchErdman");


    // Constructors

        //- Construct from dictionary
        PilchErdman(const dictionary&, CloudType&);

        //- Construct copy
        PilchErdman(const PilchErdman<CloudType>& bum);

        //- Construct and return a clone
        virtual autoPtr<BreakupModel<CloudType>> clone() const
        {
            return autoPtr<BreakupModel<CloudType>>
            (
                new PilchErdman<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~PilchErdman();


    // Member Functions

        //- Update the parcel properties
        virtual bool update
        (
            const scalar dt,
            const vector& g,
            scalar& d,
            scalar& tc,
            scalar& ms,
            scalar& nParticle,
            scalar& KHindex,
            scalar& y,
            scalar& yDot,
            const scalar d0,
            const scalar rho,
            const scalar mu,
            const scalar sigma,
            const vector& U,
            const scalar rhoc,
            const scalar muc,
            const vector& Urel,
            const scalar Urmag,
            const scalar tMom,
            scalar& dChild,
            scalar& massChild
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PilchErdman.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
