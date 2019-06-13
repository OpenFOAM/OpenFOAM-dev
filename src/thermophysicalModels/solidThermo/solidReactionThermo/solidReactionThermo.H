/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::solidReactionThermo

Description
    Foam::solidReactionThermo

SourceFiles
    solidReactionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef solidReactionThermo_H
#define solidReactionThermo_H

#include "basicSpecieMixture.H"
#include "solidThermo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class solidReactionThermo Declaration
\*---------------------------------------------------------------------------*/


class solidReactionThermo
:
    public solidThermo
{

public:

    //- Runtime type information
    TypeName("solidReactionThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidReactionThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );

     // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidReactionThermo,
        dictionary,
        (const fvMesh& mesh, const dictionary& dict, const word& phaseName),
        (mesh, dict, phaseName)
    );

    // Constructors

        //- Construct from mesh and phase name
        solidReactionThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh, dictionary and phase name
        solidReactionThermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<solidReactionThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );

         //- Standard selection based on fvMesh amd dictionary
        static autoPtr<solidReactionThermo> New
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~solidReactionThermo();


    // Member Functions

        //- Return the composition of the multi-component mixture
        virtual basicSpecieMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const basicSpecieMixture& composition() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
