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
    Foam::Function2s::binaryDiffusionCoefficient

Description
    Binary diffusion coefficient polynomial function of the log of temperature
    correcting the standard kinetic theory temperature dependence.

    \verbatim
        Dji = T^(3/2)*(c0 + c1*lnT + c2*lnT^2 + c3*lnT^3 + c4*lnT^4)/p
    \endverbatim

SourceFiles
    binaryDiffusionCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef binaryDiffusionCoefficient_H
#define binaryDiffusionCoefficient_H

#include "Function2.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{

/*---------------------------------------------------------------------------*\
          Class binaryDiffusionCoefficient Declaration
\*---------------------------------------------------------------------------*/

class binaryDiffusionCoefficient
:
    public FieldFunction2<scalar, binaryDiffusionCoefficient>
{
    // Private Data

        //- Polynomial coefficients
        Polynomial<5> coeffs_;


public:

    //- Runtime type information
    TypeName("binaryDiffusionCoefficient");


    // Constructors

        //- Construct from name and coefficients
        binaryDiffusionCoefficient
        (
            const word& name,
            const Polynomial<5> coeffs
        );

        //- Construct from name and dictionary
        binaryDiffusionCoefficient
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
                new binaryDiffusionCoefficient(*this)
            );
        }


    // Member Functions

        //- Binary diffusion coefficient function of pressure and temperature
        virtual scalar value(scalar p, scalar T) const
        {
            return T*sqrt(T)*coeffs_.value(log(T))/p;
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
