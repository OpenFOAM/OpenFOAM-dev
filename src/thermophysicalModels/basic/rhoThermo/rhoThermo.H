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
    Foam::rhoThermo

Description
    Base-class for fluid thermodynamic properties based on density.

See also
    Foam::basicThermo

SourceFiles
    rhoThermo.C

\*---------------------------------------------------------------------------*/

#ifndef rhoThermo_H
#define rhoThermo_H

#include "fluidThermo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class rhoThermo Declaration
\*---------------------------------------------------------------------------*/

class rhoThermo
:
    virtual public fluidThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;

        //- Forward declare the composite class
        class composite;


    //- Runtime type information
    TypeName("rhoThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<rhoThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~rhoThermo();


    // Member Functions

        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const = 0;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const = 0;

            //- Rename and return the thermodynamic density field [kg/m^3]
            //  This is used by solvers which create a separate continuity rho
            virtual tmp<volScalarField> renameRho() = 0;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho() = 0;

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            virtual void correctRho(const volScalarField& deltaRho) = 0;
};


/*---------------------------------------------------------------------------*\
                  Class rhoThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class rhoThermo::implementation
:
    virtual public rhoThermo
{
protected:

    // Protected data

        //- Density field [kg/m^3]
        //  Named 'rho' to avoid (potential) conflict with solver density
        volScalarField rho_;


public:

    // Constructors

        //- Construct from mesh and phase name
        implementation(const fvMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        implementation(const implementation&) = delete;


    //- Destructor
    virtual ~implementation();


    // Member Functions

        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const;

            //- Rename and return the thermodynamic density field [kg/m^3]
            //  This is used by solvers which create a separate continuity rho
            virtual tmp<volScalarField> renameRho();

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho();

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            virtual void correctRho(const volScalarField& deltaRho);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const implementation&) = delete;
};


/*---------------------------------------------------------------------------*\
                    Class rhoThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class rhoThermo::composite
:
    public basicThermo::implementation,
    public fluidThermo::implementation,
    public rhoThermo::implementation
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
            rhoThermo::implementation(mesh, phaseName)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
