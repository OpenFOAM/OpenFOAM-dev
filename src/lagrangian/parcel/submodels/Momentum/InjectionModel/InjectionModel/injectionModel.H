/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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
    Foam::injectionModel

Description
    Non-templated base class for lagrangian injection models

SourceFiles
    injectionModel.C

\*---------------------------------------------------------------------------*/

#ifndef injectionModel_H
#define injectionModel_H

#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class randomGenerator;

/*---------------------------------------------------------------------------*\
                          Class injectionModel Declaration
\*---------------------------------------------------------------------------*/

class injectionModel
{
public:

    // Public Enumerations

        //- Enumeration for the parcels' uniform size
        enum class uniformParcelSize
        {
            nParticle,
            surfaceArea,
            volume
        };

        //- Names of the parcels' uniform size
        static const NamedEnum<uniformParcelSize, 3> uniformParcelSizeNames_;


    // Static Member Functions

        //- Return a scalar uniformly distributed between zero and one. Samples
        //  on the master and then scatters the result to all processors. Use
        //  of this is discouraged. It is expensive and introduces
        //  non-randomness in all cores other then the master. It is better
        //  practice to keep a generator synchronised by calling it the same
        //  number of times on the different processes.
        static scalar globalScalar01(randomGenerator& rndGen);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
