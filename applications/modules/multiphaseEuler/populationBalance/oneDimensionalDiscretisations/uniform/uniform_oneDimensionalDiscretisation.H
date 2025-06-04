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
    Foam::oneDimensionalDiscretisations::uniform

Description
    Uniform coordinate spacing. The i-th coordinate is given by:

    \f[
        x_i^q = x_{\textrm{min}}^q + \frac{i}{n - 1} \\
        \left( x_{\textrm{max}}^q - x_{\textrm{min}}^q \right)
    \f]

Usage
    Example specification:
    \verbatim
    {
        type    uniform;

        // First and last coordinates
        min     0.01;
        max     0.51;

        // Optional moment/power of the coordinate to be spaced uniformly.
        // Defaults to one.
        q       1;
    }
    \endverbatim

SourceFiles
    oneDimensionalDiscretisation.C

\*---------------------------------------------------------------------------*/

#ifndef uniform_oneDimensionalDiscretisation_H
#define uniform_oneDimensionalDiscretisation_H

#include "oneDimensionalDiscretisation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace oneDimensionalDiscretisations
{

/*---------------------------------------------------------------------------*\
                           Class uniform Declaration
\*---------------------------------------------------------------------------*/

class uniform
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
    TypeName("uniform");


    // Public Static Member Functions

        //- Return uniform coordinates in the range 0 -> 1
        static tmp<scalarField> coordinates01(const label n);


    // Constructors

        //- Construct from a dictionary
        uniform
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
