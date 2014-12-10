/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::removeInjection

Description
    All mass available to be removed from the system is removed.

SourceFiles
    removeInjection.C

\*---------------------------------------------------------------------------*/

#ifndef removeInjection_H
#define removeInjection_H

#include "injectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class removeInjection Declaration
\*---------------------------------------------------------------------------*/

class removeInjection
:
    public injectionModel
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        removeInjection(const removeInjection&);

        //- Disallow default bitwise assignment
        void operator=(const removeInjection&);


protected:

        //- Stable film thickness - mass only removed if thickness execeeds
        //  this threhold value
        scalar deltaStable_;

        //- Mask per cell to indicate whether mass can be removed
        scalarField mask_;


public:

    //- Runtime type information
    TypeName("removeInjection");


    // Constructors

        //- Construct from surface film model
        removeInjection(surfaceFilmModel& owner, const dictionary& dict);


    //- Destructor
    virtual ~removeInjection();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct
            (
                scalarField& availableMass,
                scalarField& massToInject,
                scalarField& diameterToInject
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
