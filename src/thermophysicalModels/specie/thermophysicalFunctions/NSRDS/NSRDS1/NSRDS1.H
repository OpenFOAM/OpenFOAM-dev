/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::Function1s::NSRDS1

Description
    NSRDS function number 101

    Source:
    \verbatim
                      NSRDS - AICHE
                 Data Compilation Tables
                    of Properties of
                     Pure Compounds

        Design Institute for Physical Property Data
          American Institute of Chemical Engineers
                  345 East 47th Street
                New York, New York 10017

         National Standard Reference Data System
         American Institute of Chemical Engineers

          T.E. Daubert       -       R.P. Danner

            Department of Chemical Engineering
            The Pennsylvania State University
                University Park, PA 16802
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef NSRDS1_H
#define NSRDS1_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                           Class NSRDS1 Declaration
\*---------------------------------------------------------------------------*/

class NSRDS1
:
    public FieldFunction1<scalar, NSRDS1>
{
    // Private Data

        // NSRDS function 101 coefficients
        scalar a_, b_, c_, d_, e_;


public:

    //- Runtime type information
    TypeName("NSRDS1");


    // Constructors

        //- Construct from components
        NSRDS1
        (
            const word& name,
            const scalar a,
            const scalar b,
            const scalar c,
            const scalar d,
            const scalar e
        );

        //- Construct from name and dictionary
        NSRDS1
        (
            const word& name,
            const unitSets& units,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual tmp<Function1<scalar>> clone() const
        {
            return tmp<Function1<scalar>>(new NSRDS1(*this));
        }


    // Member Functions

        //- Return the value at a scalar coordinate
        virtual scalar value(const scalar T) const
        {
            return exp(a_ + b_/T + c_*log(T) + d_*pow(T, e_));
        }

        //- Return the derivative at a scalar coordinate
        virtual scalar derivative(const scalar T) const
        {
            return (- b_/sqr(T) + c_/T + d_*e_*pow(T, e_ - 1))*value(T);
        }

        //- Return the integral between two scalar coordinates
        virtual scalar integral(const scalar T1, const scalar T2) const
        {
            NotImplemented;
            return 0;
        }

        //- Write in dictionary format
        virtual void write(Ostream& os, const unitSets& units) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
