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
    Foam::oneDimensionalDiscretisations::exponential

Description
    Exponential coordinate spacing. The i-th coordinate is given by:

    \f[
        x_i = \left( \frac{x_{\textrm{max}}}{x_{\textrm{min}}} \\
        \right)^{\frac{i}{n - 1}} x_{\textrm{min}}
    \f]

Usage
    Example specification:
    \verbatim
    {
        type    exponential;

        // First and last coordinates
        min     0.01;
        max     0.51;
    }
    \endverbatim

SourceFiles
    oneDimensionalDiscretisation.C

\*---------------------------------------------------------------------------*/

#ifndef exponential_oneDimensionalDiscretisation_H
#define exponential_oneDimensionalDiscretisation_H

#include "oneDimensionalDiscretisation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace oneDimensionalDiscretisations
{

/*---------------------------------------------------------------------------*\
                          Class exponential Declaration
\*---------------------------------------------------------------------------*/

class exponential
:
    public oneDimensionalDiscretisation
{
private:

    // Private Static Member Functions

        //- Generate the coordinates
        tmp<scalarField> coordinates
        (
            const dimensionSet& dims,
            const label n,
            const dictionary& dict
        );


public:

    //- Runtime type information
    TypeName("exponential");


    // Constructors

        //- Construct from a dictionary
        exponential
        (
            const word& name,
            const dimensionSet& dims,
            const label n,
            const dictionary& dict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace oneDimensionalDiscretisations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
