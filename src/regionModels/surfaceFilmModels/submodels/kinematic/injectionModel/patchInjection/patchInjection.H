/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::patchInjection

Description
    Remove and inject the mass in the film as it passes over the selected
    patches.

SourceFiles
    patchInjection.C

\*---------------------------------------------------------------------------*/

#ifndef patchInjection_H
#define patchInjection_H

#include "injectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class patchInjection Declaration
\*---------------------------------------------------------------------------*/

class patchInjection
:
    public injectionModel
{
protected:

        //- Stable film thickness - mass only removed if thickness exceeds
        //  this threshold value
        scalar deltaStable_;

        //- List of patch IDs at which the film is removed
        labelList patchIDs_;

        //- Injected mass for each patch at which the film is removed
        scalarField patchInjectedMasses_;


public:

    //- Runtime type information
    TypeName("patchInjection");


    // Constructors

        //- Construct from surface film model
        patchInjection(surfaceFilmRegionModel& film, const dictionary& dict);

        //- Disallow default bitwise copy construction
        patchInjection(const patchInjection&) = delete;


    //- Destructor
    virtual ~patchInjection();


    // Member Functions

        //- Correct
        virtual void correct
        (
            scalarField& availableMass,
            scalarField& massToInject,
            scalarField& diameterToInject
        );

        //- Accumulate the total mass injected for the patches into the
        //  scalarField provided
        virtual void patchInjectedMassTotals
        (
            scalarField& patchMasses
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const patchInjection&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
