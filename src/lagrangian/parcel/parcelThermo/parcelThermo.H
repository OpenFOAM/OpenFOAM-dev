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
    Foam::parcelThermo

Description
    Thermo package for (S)olids (L)iquids and (G)ases
    Takes reference to thermo package, and provides:
    - liquids : liquid components - access  to elemental properties
    - solids  : solid components - access  to elemental properties

    If no liquids or solids are specified, their respective pointers will also
    be nullptr.

    Registered to the mesh so that it can be looked-up

SourceFiles
    parcelThermo.C

\*---------------------------------------------------------------------------*/

#ifndef parcelThermo_H
#define parcelThermo_H

#include "fluidThermo.H"
#include "liquidMixtureProperties.H"
#include "solidMixtureProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class parcelThermo Declaration
\*---------------------------------------------------------------------------*/

class parcelThermo
{
    // Private Data

        //- Additional liquid properties data
        autoPtr<liquidMixtureProperties> liquids_;

        //- Additional solid properties data
        autoPtr<solidMixtureProperties> solids_;


public:

    // Constructors

        //- Construct from carrier thermo
        parcelThermo(const fluidThermo& carrierThermo);


    //- Destructor
    virtual ~parcelThermo();


    // Member Functions

        //- Return reference to the global (additional) liquids
        const liquidMixtureProperties& liquids() const;

        //- Return reference to the global (additional) solids
        const solidMixtureProperties& solids() const;

        //- Index of liquid component
        label liquidId
        (
            const word& cmptName,
            bool allowNotFound = false
        ) const;

        //- Index of solid component
        label solidId
        (
            const word& cmptName,
            bool allowNotFound = false
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
