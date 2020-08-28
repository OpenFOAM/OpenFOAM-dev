/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::TAB

Description
    The TAB Method for Numerical Calculation of Spray Droplet Breakup.

    @verbatim
        O'Rourke, P.J. and Amsden, A.A.,
        "The TAB Method for Numerical Calculation of Spray Droplet Breakup,"
        1987 SAE International Fuels and Lubricants Meeting and Exposition,
        Toronto, Ontario, November 2-5, 1987,
        Los Alamos National Laboratory document LA-UR-87-2105;
        SAE Technical Paper Series, Paper 872089.
    @endverbatim

    This implementation follows the kiva version.

See also
    The Enhanced %TAB model - ETAB

\*---------------------------------------------------------------------------*/

#ifndef TAB_H
#define TAB_H

#include "BreakupModel.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                             Class TAB Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TAB
:
    public BreakupModel<CloudType>
{
public:

    //- Enumeration for the SMD breakup calculation
    enum SMDMethods
    {
        method1,
        method2
    };


private:

    // Private Data

        // Inverse function approximation of the Rossin-Rammler Distribution
        // used when calculating the droplet size after breakup
        FixedList<scalar, 100> rrd_;


        // Model constants

            word SMDCalcMethod_;
            SMDMethods SMDMethod_;


public:

    //- Runtime type information
    TypeName("TAB");


    // Constructors

        //- Construct from dictionary
        TAB(const dictionary& dict, CloudType& owner);

        //- Construct copy
        TAB(const TAB<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<BreakupModel<CloudType>> clone() const
        {
            return autoPtr<BreakupModel<CloudType>>
            (
                new TAB<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~TAB();


    // Member Functions

        //- Update the parcel diameter
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
    #include "TAB.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
