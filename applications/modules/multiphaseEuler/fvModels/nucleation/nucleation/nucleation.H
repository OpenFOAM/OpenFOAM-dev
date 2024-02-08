/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::fv::nucleation

Description
    Mix-in interface for nucleation models. Provides access to properties of
    the nucleation process, such as diameter and rate of production of nuclei.

SourceFiles
    nucleation.C

\*---------------------------------------------------------------------------*/

#ifndef nucleation_H
#define nucleation_H

#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class GeoMesh>
class DimensionedField;

namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class nucleation Declaration
\*---------------------------------------------------------------------------*/

class nucleation
{
public:

    //- Runtime type information
    TypeName("nucleation");


    // Constructors

        //- Default construct
        nucleation();


    // Destructor
    virtual ~nucleation();


    // Member Functions

        // Sources

            //- Return the diameter of nuclei
            virtual tmp<DimensionedField<scalar, volMesh>> d() const = 0;

            //- Return the number rate at which nuclei are generated
            virtual tmp<DimensionedField<scalar, volMesh>> nDot() const = 0;

            //- Return the mass transfer rate
            virtual tmp<DimensionedField<scalar, volMesh>> mDot() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
