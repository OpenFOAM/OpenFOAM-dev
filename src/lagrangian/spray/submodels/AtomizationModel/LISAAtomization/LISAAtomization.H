/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    Foam::LISAAtomization

Description
    Primary Breakup Model for pressure swirl atomizers.

    Accurate description in
    @verbatim
    P.K. Senecal, D.P. Schmidt, I. Nouar, C.J. Rutland, R.D. Reitz, M. Corradini
    "Modeling high-speed viscous liquid sheet atomization"
    International Journal of Multiphase Flow 25 (1999) pags. 1073-1097
    @endverbatim

    and

    @verbatim
    D.P. Schmidt, I. Nouar, P.K. Senecal, C.J. Rutland, J.K. Martin, R.D. Reitz
    "Pressure-Swirl Atomization in the Near Field"
    SAE Techical Paper Series 1999-01-0496
    @endverbatim

\*---------------------------------------------------------------------------*/

#ifndef LISAAtomization_H
#define LISAAtomization_H

#include "AtomizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class LISAAtomization Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class LISAAtomization
:
    public AtomizationModel<CloudType>
{

public:

    //- Enumeration for SMD calculations
    enum SMDMethods
    {
        method1,
        method2
    };


private:

    // Private data

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
        LISAAtomization(const dictionary&, CloudType&);

        //- Construct copy
        LISAAtomization(const LISAAtomization<CloudType>& am);

        //- Construct and return a clone
        virtual autoPtr<AtomizationModel<CloudType>> clone() const
        {
            return autoPtr<AtomizationModel<CloudType>>
            (
                new LISAAtomization<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~LISAAtomization();


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
    #include "LISAAtomization.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
