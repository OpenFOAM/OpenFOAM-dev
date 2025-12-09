/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::Function2s::FullerDiffusionCoefficient

Description
    Fuller et al. correlation for vapour mass diffusivity.

    \verbatim
        E. N. Fuller, E. N., Schettler P. D.,  Giddings, J. C. (1966).
        A New Method For Prediction Of Binary Gas-Phase Diffusion Coefficients
        Ind. Eng. Chem. 58(5), 19.
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef FullerDiffusionCoefficient_H
#define FullerDiffusionCoefficient_H

#include "Function2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{

/*---------------------------------------------------------------------------*\
                         Class FullerDiffusionCoefficient Declaration
\*---------------------------------------------------------------------------*/

class FullerDiffusionCoefficient
:
    public FieldFunction2<scalar, FullerDiffusionCoefficient>
{
    // Private Data

        // Vapour mass diffusivity function coefficients
        scalar Va_, Vb_, Ma_, Mb_;

        // Helper variables
        scalar alpha_, beta_;


public:

    //- Runtime type information
    TypeName("FullerDiffusionCoefficient");


    // Constructors

        //- Construct from name and dictionary
        FullerDiffusionCoefficient
        (
            const word& name,
            const unitConversions& units,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual tmp<Function2<scalar>> clone() const
        {
            return tmp<Function2<scalar>>
            (
                new FullerDiffusionCoefficient(*this)
            );
        }


    // Member Functions

        //- Inherit base class value methods to prevent clang warnings
        using FieldFunction2<scalar, FullerDiffusionCoefficient>::value;

        //- Vapour mass diffusivity function
        virtual scalar value(scalar p, scalar T) const
        {
            return 101325e-7*(pow(T, 1.75))*alpha_/(p*beta_);
        }

        //- Write the function coefficients
        virtual void write(Ostream& os, const unitConversions& units) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function2s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
