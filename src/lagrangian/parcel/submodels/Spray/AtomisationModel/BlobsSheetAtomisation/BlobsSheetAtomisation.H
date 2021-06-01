/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::BlobsSheetAtomisation

Description
    Primary Breakup Model for pressure swirl atomisers.

    Accurate description in
    @verbatim
    Z. Han, S. Parrish, P.V. Farrell, R.D. Reitz
    "Modeling Atomisation Processes Of Pressure Swirl Hollow-Cone Fuel Sprays"
    Atomisation and Sprays, vol. 7, pp. 663-684, 1997

    and

    L. Allocca, G. Bella, A. De Vita, L. Di Angelo
    "Experimental Validation of a GDI Spray Model"
    SAE Technical Paper Series, 2002-01-1137
    @endverbatim


\*---------------------------------------------------------------------------*/

#ifndef BlobsSheetAtomisation_H
#define BlobsSheetAtomisation_H

#include "AtomisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class BlobsSheetAtomisation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BlobsSheetAtomisation
:
    public AtomisationModel<CloudType>
{
    scalar B_;
    scalar angle_;


public:

    //- Runtime type information
    TypeName("blobsSheetAtomisation");


    // Constructors

        //- Construct from dictionary
        BlobsSheetAtomisation(const dictionary& dict, CloudType& cloud);

        //- Construct copy
        BlobsSheetAtomisation(const BlobsSheetAtomisation<CloudType>& am);

        //- Construct and return a clone
        virtual autoPtr<AtomisationModel<CloudType>> clone() const
        {
            return autoPtr<AtomisationModel<CloudType>>
            (
                new BlobsSheetAtomisation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~BlobsSheetAtomisation();


    // Member Functions

        //- Initial value of liquidCore
        virtual scalar initLiquidCore() const;

        //- Flag to indicate if chi needs to be calculated
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
            Random& rndGen
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BlobsSheetAtomisation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
