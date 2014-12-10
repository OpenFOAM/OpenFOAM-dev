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
    Foam::BlobsSheetAtomization

Description
    Primary Breakup Model for pressure swirl atomizers.

    Accurate description in
    @verbatim
    Z. Han, S. Parrish, P.V. Farrell, R.D. Reitz
    "Modeling Atomization Processes Of Pressure Swirl Hollow-Cone Fuel Sprays"
    Atomization and Sprays, vol. 7, pp. 663-684, 1997

    and

    L. Allocca, G. Bella, A. De Vita, L. Di Angelo
    "Experimental Validation of a GDI Spray Model"
    SAE Technical Paper Series, 2002-01-1137
    @endverbatim


\*---------------------------------------------------------------------------*/

#ifndef BlobsSheetAtomization_H
#define BlobsSheetAtomization_H

#include "AtomizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                   Class BlobsSheetAtomization Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BlobsSheetAtomization
:
    public AtomizationModel<CloudType>
{
private:

    scalar B_;
    scalar angle_;


public:

    //- Runtime type information
    TypeName("blobsSheetAtomization");


    // Constructors

        //- Construct from dictionary
        BlobsSheetAtomization(const dictionary& dict, CloudType& cloud);

        //- Construct copy
        BlobsSheetAtomization(const BlobsSheetAtomization<CloudType>& am);

        //- Construct and return a clone
        virtual autoPtr<AtomizationModel<CloudType> > clone() const
        {
            return autoPtr<AtomizationModel<CloudType> >
            (
                new BlobsSheetAtomization<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~BlobsSheetAtomization();


    // Member Functions

        //- initial value of liquidCore
        virtual scalar initLiquidCore() const;

        //- flag to indicate if chi needs to be calculated
        virtual bool calcChi() const;

        virtual void update
        (
            const scalar dt,
            scalar& d,
            scalar& liquidCore,
            scalar& tc,
            const scalar rho,
            const scalar mu,
            const scalar sigma,
            const scalar volFlowRate,
            const scalar rhoAv,
            const scalar Urel,
            const vector& pos,
            const vector& injectionPos,
            const scalar pAmbient,
            const scalar chi,
            cachedRandom& rndGen
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlobsSheetAtomization.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
