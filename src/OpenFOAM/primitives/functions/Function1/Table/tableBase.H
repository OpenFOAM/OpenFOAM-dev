/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::tableBase

SourceFiles
    tableBase.C

\*---------------------------------------------------------------------------*/

#ifndef tableBase_H
#define tableBase_H

#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                          Class tableBase Declaration
\*---------------------------------------------------------------------------*/

class tableBase
{
public:

    // Public data types

        //- Enumeration for handling out-of-bound values
        enum class boundsHandling
        {
            error,  // Exit with a FatalError
            warn,   // Issue warning and clamp value (default)
            clamp,  // Clamp value to the start/end value
            repeat  // Treat as a repeating list
        };

        //- Enumeration names for handling out-of-bound values
        static const NamedEnum<boundsHandling, 4> boundsHandlingNames_;


    // Constructors

        //- Construct null
        tableBase()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
