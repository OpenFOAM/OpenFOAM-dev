/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::psiReactionThermo

Description
    Base-class for multi-component fluid thermodynamic properties based on
    compressibility.

See also
    Foam::basicThermo

SourceFiles
    psiReactionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef psiReactionThermo_H
#define psiReactionThermo_H

#include "psiThermo.H"
#include "fluidReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class psiReactionThermo Declaration
\*---------------------------------------------------------------------------*/

class psiReactionThermo
:
    virtual public psiThermo,
    virtual public fluidReactionThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;

        //- Forward declare the composite class
        class composite;


    //- Runtime type information
    TypeName("psiReactionThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        psiReactionThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<psiReactionThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~psiReactionThermo();
};


/*---------------------------------------------------------------------------*\
              Class psiReactionThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class psiReactionThermo::implementation
:
    virtual public psiReactionThermo
{
public:

    // Constructors

        //- Construct from mesh and phase name
        implementation(const fvMesh&, const word& phaseName);


    //- Destructor
    virtual ~implementation();
};


/*---------------------------------------------------------------------------*\
                 Class psiReactionThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class psiReactionThermo::composite
:
    public basicThermo::implementation,
    public fluidThermo::implementation,
    public psiThermo::implementation,
    public fluidReactionThermo::implementation,
    public psiReactionThermo::implementation
{
public:

    // Constructors

        //- Construct from mesh and phase name
        composite
        (
            const fvMesh& mesh,
            const word& phaseName
        )
        :
            basicThermo::implementation(mesh, phaseName),
            fluidThermo::implementation(mesh, phaseName),
            psiThermo::implementation(mesh, phaseName),
            fluidReactionThermo::implementation(mesh, phaseName),
            psiReactionThermo::implementation(mesh, phaseName)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
