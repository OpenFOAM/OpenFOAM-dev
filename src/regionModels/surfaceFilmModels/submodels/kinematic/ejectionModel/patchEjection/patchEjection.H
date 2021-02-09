/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::patchEjection

Description
    Remove and eject the mass in the film as it passes over the selected
    patches.

SourceFiles
    patchEjection.C

\*---------------------------------------------------------------------------*/

#ifndef patchEjection_H
#define patchEjection_H

#include "ejectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                       Class patchEjection Declaration
\*---------------------------------------------------------------------------*/

class patchEjection
:
    public ejectionModel
{
protected:

        //- Stable film thickness - mass only removed if thickness exceeds
        //  this threshold value
        scalar deltaStable_;

        //- List of patch IDs at which the film is removed
        labelList patchIDs_;

        //- Ejected mass for each patch at which the film is removed
        scalarField patchEjectedMasses_;


public:

    //- Runtime type information
    TypeName("patchEjection");


    // Constructors

        //- Construct from surface film model
        patchEjection(surfaceFilmRegionModel& film, const dictionary& dict);

        //- Disallow default bitwise copy construction
        patchEjection(const patchEjection&) = delete;


    //- Destructor
    virtual ~patchEjection();


    // Member Functions

        //- Correct
        virtual void correct
        (
            scalarField& availableMass,
            scalarField& massToEject,
            scalarField& diameterToEject
        );

        //- Accumulate the total mass ejected for the patches into the
        //  scalarField provided
        virtual void patchEjectedMassTotals
        (
            scalarField& patchMasses
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const patchEjection&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
