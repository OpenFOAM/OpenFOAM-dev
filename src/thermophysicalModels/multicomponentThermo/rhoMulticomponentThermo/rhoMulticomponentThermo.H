/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::rhoMulticomponentThermo

Description
    Base-class for multi-component fluid thermodynamic properties based on
    density.

See also
    Foam::basicThermo

SourceFiles
    rhoMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef rhoMulticomponentThermo_H
#define rhoMulticomponentThermo_H

#include "rhoThermo.H"
#include "fluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class rhoMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

class rhoMulticomponentThermo
:
    virtual public rhoThermo,
    virtual public fluidMulticomponentThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;

        //- Forward declare the composite class
        class composite;


    //- Runtime type information
    TypeName("rhoMulticomponentThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoMulticomponentThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<rhoMulticomponentThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~rhoMulticomponentThermo();
};


/*---------------------------------------------------------------------------*\
              Class rhoMulticomponentThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class rhoMulticomponentThermo::implementation
:
    virtual public rhoMulticomponentThermo
{

public:

    // Constructors

        //- Construct from mesh and phase name
        implementation(const fvMesh&, const word& phaseName);


    //- Destructor
    virtual ~implementation();
};


/*---------------------------------------------------------------------------*\
                Class rhoMulticomponentThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class rhoMulticomponentThermo::composite
:
    public basicThermo::implementation,
    public fluidThermo::implementation,
    public rhoThermo::implementation,
    public fluidMulticomponentThermo::implementation,
    public rhoMulticomponentThermo::implementation
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
            rhoThermo::implementation(mesh, phaseName),
            fluidMulticomponentThermo::implementation(mesh, phaseName),
            rhoMulticomponentThermo::implementation(mesh, phaseName)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
