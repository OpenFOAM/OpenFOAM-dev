/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::enginePiston

Description
    Foam::enginePiston

SourceFiles
    enginePiston.C

\*---------------------------------------------------------------------------*/

#ifndef enginePiston_H
#define enginePiston_H

#include "fvMeshMoversEngine.H"
#include "polyPatchDynamicID.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class enginePiston Declaration
\*---------------------------------------------------------------------------*/

class enginePiston
{
    // Private Data

        //- Reference to engine mesh mover
        const fvMeshMovers::engine& meshMover_;

        //- Piston patch
        polyPatchDynamicID patchID_;

        //- Coordinate system
        autoPtr<coordinateSystem> csPtr_;


        // Piston layering data

            //- Min layer thickness
            const scalar minLayer_;

            //- Max layer thickness
            const scalar maxLayer_;


public:

    // Constructors

        //- Construct from components
        enginePiston
        (
            const fvMeshMover& meshMover,
            const word& pistonPatchName,
            const autoPtr<coordinateSystem>& pistonCS,
            const scalar minLayer,
            const scalar maxLayer

        );

        //- Construct from dictionary
        enginePiston
        (
            const fvMeshMover& meshMover,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        enginePiston(const enginePiston&) = delete;


    // Member Functions

        //- Return coordinate system
        const coordinateSystem& cs() const
        {
            return csPtr_();
        }

        //- Return ID of piston patch
        const polyPatchDynamicID& patchID() const
        {
            return patchID_;
        }

        // Piston layering thickness

            scalar minLayer() const
            {
                return minLayer_;
            }

            scalar maxLayer() const
            {
                return maxLayer_;
            }


        //- Write dictionary
        void writeDict(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const enginePiston&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
