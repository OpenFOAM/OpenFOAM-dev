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
    Foam::Function2s::APIdiffCoefFunc

Description
    API function for vapour mass diffusivity

    Source:
    \verbatim
            API (American Petroleum Institute)
                    Technical Data Book
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef APIdiffCoef_H
#define APIdiffCoef_H

#include "Function2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{

/*---------------------------------------------------------------------------*\
                         Class APIdiffCoef Declaration
\*---------------------------------------------------------------------------*/

class APIdiffCoef
:
    public FieldFunction2<scalar, APIdiffCoef>
{
    // Private Data

        // API vapour mass diffusivity function coefficients
        scalar a_, b_, wf_, wa_;

        // Helper variables
        scalar alpha_, beta_;


public:

    //- Runtime type information
    TypeName("APIdiffCoef");


    // Constructors

        //- Construct from components
        APIdiffCoef
        (
            const word& name,
            const scalar a,
            const scalar b,
            const scalar wf,
            const scalar wa
        );

        //- Construct from name and dictionary
        APIdiffCoef
        (
            const word& name,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual tmp<Function2<scalar>> clone() const
        {
            return tmp<Function2<scalar>>(new APIdiffCoef(*this));
        }


    // Member Functions

        //- Inherit base class value methods to prevent clang warnings
        using FieldFunction2<scalar, APIdiffCoef>::value;

        //- API vapour mass diffusivity function using properties from
        //  construction
        virtual scalar value(scalar p, scalar T) const
        {
            return 3.6059e-3*(pow(1.8*T, 1.75))*alpha_/(p*beta_);
        }

        //- API vapour mass diffusivity function using properties from
        //  construction - with specified binary pair
        scalar value(scalar p, scalar T, scalar  Wa) const
        {
            const scalar alphaBinary = sqrt(1/wf_ + 1/Wa);
            return 3.6059e-3*(pow(1.8*T, 1.75))*alphaBinary/(p*beta_);
        }

        //- Write the function coefficients
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function2s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
