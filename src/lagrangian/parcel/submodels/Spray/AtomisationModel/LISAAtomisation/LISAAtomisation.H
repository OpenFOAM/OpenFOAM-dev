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
    Foam::LISAAtomisation

Description
    Primary Breakup Model for pressure swirl atomisers.

    Accurate description in
    @verbatim
    P.K. Senecal, D.P. Schmidt, I. Nouar, C.J. Rutland, R.D. Reitz, M. Corradini
    "Modeling high-speed viscous liquid sheet atomisation"
    International Journal of Multiphase Flow 25 (1999) pags. 1073-1097
    @endverbatim

    and

    @verbatim
    D.P. Schmidt, I. Nouar, P.K. Senecal, C.J. Rutland, J.K. Martin, R.D. Reitz
    "Pressure-Swirl Atomisation in the Near Field"
    SAE Techical Paper Series 1999-01-0496
    @endverbatim

\*---------------------------------------------------------------------------*/

#ifndef LISAAtomisation_H
#define LISAAtomisation_H

#include "AtomisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class LISAAtomisation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class LISAAtomisation
:
    public AtomisationModel<CloudType>
{

public:

    //- Enumeration for SMD calculations
    enum SMDMethods
    {
        method1,
        method2
    };


private:

    // Private Data

        scalar Cl_;
        scalar cTau_;
        scalar Q_;
        scalar lisaExp_;
        vector injectorDirection_;
        word SMDCalcMethod_;

        SMDMethods SMDMethod_;


public:

    //- Runtime type information
    TypeName("LISA");


    // Constructors

        //- Construct from dictionary
        LISAAtomisation(const dictionary&, CloudType&);

        //- Construct copy
        LISAAtomisation(const LISAAtomisation<CloudType>& am);

        //- Construct and return a clone
        virtual autoPtr<AtomisationModel<CloudType>> clone() const
        {
            return autoPtr<AtomisationModel<CloudType>>
            (
                new LISAAtomisation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~LISAAtomisation();


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
    #include "LISAAtomisation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
