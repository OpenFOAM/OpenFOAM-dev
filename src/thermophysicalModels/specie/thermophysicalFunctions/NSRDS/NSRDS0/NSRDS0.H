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
    Foam::Function1s::NSRDS0

Description
    NSRDS function number 100

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

#ifndef NSRDS0_H
#define NSRDS0_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                           Class NSRDS0 Declaration
\*---------------------------------------------------------------------------*/

class NSRDS0
:
    public FieldFunction1<scalar, NSRDS0>
{
    // Private Data

        // NSRDS function 100 coefficients
        scalar a_, b_, c_, d_, e_, f_;


public:

    //- Runtime type information
    TypeName("NSRDS0");


    // Constructors

        //- Construct from components
        NSRDS0
        (
            const word& name,
            const scalar a,
            const scalar b,
            const scalar c,
            const scalar d,
            const scalar e,
            const scalar f
        );

        //- Construct from name and dictionary
        NSRDS0
        (
            const word& name,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual tmp<Function1<scalar>> clone() const
        {
            return tmp<Function1<scalar>>(new NSRDS0(*this));
        }


    // Member Functions

        //- Evaluate the function and return the result
        virtual scalar value(scalar T) const
        {
            return ((((f_*T + e_)*T + d_)*T + c_)*T + b_)*T + a_;
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
