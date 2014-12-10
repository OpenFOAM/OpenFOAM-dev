/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::incompressible::LESdelta

Description
    Abstract base class for incompressible LES deltas

SourceFiles
    incompressibleLESdelta.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleLESdelta_H
#define incompressibleLESdelta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// To avoid macro problems typedef scoped class
typedef Foam::LESdelta foamLESdelta;

namespace incompressible
{

/*---------------------------------------------------------------------------*\
                          Class LESdelta Declaration
\*---------------------------------------------------------------------------*/

class LESdelta
:
    public foamLESdelta
{

public:

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            foamLESdelta,
            dictionary,
            (
                const word& name,
                const fvMesh& mesh,
                const dictionary& LESdeltaDict
            ),
            (name, mesh, LESdeltaDict)
        );

    // Constructors

        //- Construct from name and mesh
        LESdelta(const word& name, const fvMesh& mesh);


    // Selectors

        //- Return a reference to the selected LES delta
        static autoPtr<foamLESdelta> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LESdelta()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
