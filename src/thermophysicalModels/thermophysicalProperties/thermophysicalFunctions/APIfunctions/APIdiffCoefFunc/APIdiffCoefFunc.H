/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Foam::APIdiffCoefFunc

Description
    API function for vapour mass diffusivity

    Source:
    \verbatim
            API (American Petroleum Institute)
                    Technical Data Book
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef APIdiffCoefFunc_H
#define APIdiffCoefFunc_H

#include "thermophysicalFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class APIdiffCoefFunc Declaration
\*---------------------------------------------------------------------------*/

class APIdiffCoefFunc
:
    public thermophysicalFunction
{
    // Private data

        // API vapour mass diffusivity function coefficients
        scalar a_, b_, wf_, wa_;

        // Helper variables
        scalar alpha_, beta_;


public:

    //- Runtime type information
    TypeName("APIdiffCoefFunc");


    // Constructors

        //- Construct from components
        APIdiffCoefFunc
        (
            const scalar a,
            const scalar b,
            const scalar wf,
            const scalar wa
        );

        //- Construct from dictionary
        APIdiffCoefFunc(const dictionary& dict);


    // Member Functions

        //- API vapour mass diffusivity function using properties from
        //  construction
        scalar f(scalar p, scalar T) const
        {
            return 3.6059e-3*(pow(1.8*T, 1.75))*alpha_/(p*beta_);
        }

        //- API vapour mass diffusivity function using properties from
        //  construction - with specified binary pair
        scalar f(scalar p, scalar T, scalar  Wa) const
        {
            const scalar alphaBinary = sqrt(1/wf_ + 1/Wa);
            return 3.6059e-3*(pow(1.8*T, 1.75))*alphaBinary/(p*beta_);
        }

        //- Write the function coefficients
        void writeData(Ostream& os) const
        {
            os  << a_ << token::SPACE
                << b_ << token::SPACE
                << wf_ << token::SPACE
                << wa_;
        }


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const APIdiffCoefFunc& f)
        {
            f.writeData(os);
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
