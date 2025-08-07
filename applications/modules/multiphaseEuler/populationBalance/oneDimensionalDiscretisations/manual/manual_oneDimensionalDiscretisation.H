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
    Foam::oneDimensionalDiscretisations::manual

Description
    Manual coordinate spacing. The user explicitly provides a list of
    coordinates which must increase monotonically and be of the required
    length.

Usage
    Example specification:
    \verbatim
    {
        type    manual;

        // Coordinates
        values  (0.01 0.11 0.21 0.31 0.41 0.51);
    }
    \endverbatim

SourceFiles
    oneDimensionalDiscretisation.C

\*---------------------------------------------------------------------------*/

#ifndef manual_oneDimensionalDiscretisation_H
#define manual_oneDimensionalDiscretisation_H

#include "oneDimensionalDiscretisation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace oneDimensionalDiscretisations
{

/*---------------------------------------------------------------------------*\
                          Class manual Declaration
\*---------------------------------------------------------------------------*/

class manual
:
    public oneDimensionalDiscretisation
{
public:

    //- Runtime type information
    TypeName("manual");


    // Constructors

        //- Construct from a dictionary
        manual
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
