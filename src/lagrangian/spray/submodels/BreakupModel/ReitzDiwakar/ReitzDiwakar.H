/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::ReitzDiwakar

Description
    secondary breakup model

    @verbatim
    Reitz, R.D.
    "Modelling atomization processes in highpressure vaporizing sprays"
    Atomization and Spray Technology 3 (1987), 309-337
    @endverbatim

    @verbatim
    Reitz, R.D. and Diwakar, R.
    "Effect of drop breakup on fuel sprays"
    SAE Tech. paper series, 860469 (1986)
    @endverbatim

    @verbatim
    Reitz, R.D. and Diwakar, R.
    "Structure of high-pressure fuel sprays"
    SAE Tech. paper series, 870598 (1987)
    @endverbatim

\*---------------------------------------------------------------------------*/

#ifndef ReitzDiwakar_H
#define ReitzDiwakar_H

#include "BreakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class ReitzDiwakar Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ReitzDiwakar
:
    public BreakupModel<CloudType>
{
private:

    // Private data

        scalar Cbag_;
        scalar Cb_;
        scalar Cstrip_;
        scalar Cs_;


public:

    //- Runtime type information
    TypeName("ReitzDiwakar");


    // Constructors

        //- Construct from dictionary
        ReitzDiwakar(const dictionary&, CloudType&);

        //- Construct copy
        ReitzDiwakar(const ReitzDiwakar<CloudType>& bum);

        //- Construct and return a clone
        virtual autoPtr<BreakupModel<CloudType> > clone() const
        {
            return autoPtr<BreakupModel<CloudType> >
            (
                new ReitzDiwakar<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ReitzDiwakar();


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
#   include "ReitzDiwakar.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
