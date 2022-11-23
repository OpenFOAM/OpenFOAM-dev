/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::surfaceFilmModels::speciePhaseChange

Description
    Specie phase change model abstract base class

SourceFiles
    speciePhaseChange.C

\*---------------------------------------------------------------------------*/

#ifndef speciePhaseChange_H
#define speciePhaseChange_H

#include "phaseChangeModel.H"
#include "thermoSurfaceFilm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                     Class speciePhaseChange Declaration
\*---------------------------------------------------------------------------*/

class speciePhaseChange
:
    public phaseChangeModel
{

public:

    //- Runtime type information
    TypeName("speciePhaseChange");


    // Constructors

        //- Construct from surface film model
        speciePhaseChange
        (
            const word& modelType,
            surfaceFilm& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        speciePhaseChange(const speciePhaseChange&) = delete;


    //- Destructor
    virtual ~speciePhaseChange();


    // Member Functions

        label vapId() const;

        scalar Wvap() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const speciePhaseChange&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
