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
    Foam::ETAB

Description
    The Enhanced TAB model.

    Described in the papers below.
    @verbatim
    F.X. Tanner
        "Liquid Jet Atomization and Droplet Breakup Modeling of
        Non-Evaporating Diesel Fuel Sprays"
        SAE 970050,
        SAE Transactions: Journal of Engines, Vol 106, Sec 3 pp 127-140

    F.X. Tanner and G. Weisser
        "Simulation of Liquid Jet Atomization for
        Fuel Sprays by Means of Cascade Drop Breakup Model"
        SAE 980808
        SAE Technical Paper Series
    @endverbatim

See Also
    The TAB model

\*---------------------------------------------------------------------------*/

#ifndef ETAB_H
#define ETAB_H

#include "BreakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class ETAB Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ETAB
:
    public BreakupModel<CloudType>
{
private:

    // Private data

        // Model constants

            scalar k1_;
            scalar k2_;
            scalar WeTransition_;
            scalar AWe_;


public:

    //- Runtime type information
    TypeName("ETAB");


    // Constructors

        //- Construct from dictionary
        ETAB(const dictionary&, CloudType&);

        //- Construct copy
        ETAB(const ETAB<CloudType>& bum);

        //- Construct and return a clone
        virtual autoPtr<BreakupModel<CloudType> > clone() const
        {
            return autoPtr<BreakupModel<CloudType> >
            (
                new ETAB<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ETAB();


    // Member Functions

        //- update the parcel properties
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
#   include "ETAB.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
