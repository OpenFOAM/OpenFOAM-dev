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
    Foam::Function1s::NSRDS7

Description
    NSRDS-AICHE function number 107

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

#ifndef NSRDS7_H
#define NSRDS7_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                           Class NSRDS7 Declaration
\*---------------------------------------------------------------------------*/

class NSRDS7
:
    public FieldFunction1<scalar, NSRDS7>
{
    // Private Data

        // NSRDS function 107 coefficients
        scalar a_, b_, c_, d_, e_;


public:

    //- Runtime type information
    TypeName("NSRDS7");


    // Constructors

        //- Construct from components
        NSRDS7
        (
            const word& name,
            const scalar a,
            const scalar b,
            const scalar c,
            const scalar d,
            const scalar e
        );

        //- Construct from name and dictionary
        NSRDS7
        (
            const word& name,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual tmp<Function1<scalar>> clone() const
        {
            return tmp<Function1<scalar>>(new NSRDS7(*this));
        }


    // Member Functions

        //- Evaluate the function and return the result
        virtual scalar value(scalar T) const
        {
            return a_ + b_*sqr((c_/T)/sinh(c_/T)) + d_*sqr((e_/T)/cosh(e_/T));
        }

        //- Integrate between two scalar values
        virtual scalar integral(const scalar x1, const scalar x2) const;

        //- Write the function coefficients
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
