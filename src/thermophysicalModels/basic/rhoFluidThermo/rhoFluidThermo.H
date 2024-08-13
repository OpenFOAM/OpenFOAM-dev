/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::rhoFluidThermo

Description
    Base-class for fluid thermodynamic properties based on density.

See also
    Foam::basicThermo

SourceFiles
    rhoFluidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef rhoFluidThermo_H
#define rhoFluidThermo_H

#include "RhoFluidThermo.H"
#include "pureThermo.H"
#include "rhoThermo.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class rhoFluidThermo Declaration
\*---------------------------------------------------------------------------*/

class rhoFluidThermo
:
    virtual public rhoThermo,
    virtual public fluidThermo
{
public:

    // Public Classes

        //- Forward declare the composite class
        class composite;


    // Public Typedefs

        //- The derived type
        template<class MixtureType>
        using DerivedThermoType =
            RhoFluidThermo<BasicThermo<MixtureType, composite>>;

        //- The derived name
        static const word derivedThermoName;


    //- Runtime type information
    TypeName("rhoFluidThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoFluidThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<rhoFluidThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~rhoFluidThermo();


    // Member Functions

        // Derived Thermodynamic Properties

            //- Rename and return the thermodynamic density field [kg/m^3]
            //  This is used by solvers which create a separate continuity rho
            virtual tmp<volScalarField> renameRho();

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            virtual void correctRho(const volScalarField& deltaRho);
};


/*---------------------------------------------------------------------------*\
                    Class rhoFluidThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class rhoFluidThermo::composite
:
    public basicThermo::implementation,
    public pureThermo,
    public rhoThermo::implementation,
    public fluidThermo::implementation,
    public rhoFluidThermo
{
public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        composite
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& phaseName
        )
        :
            basicThermo::implementation(dict, mesh, phaseName),
            rhoThermo::implementation(dict, mesh, phaseName),
            fluidThermo::implementation(dict, mesh, phaseName)
        {}

        //- Construct from dictionary, mesh and phase name
        template<class MixtureType>
        composite
        (
            const dictionary& dict,
            const MixtureType& mixture,
            const fvMesh& mesh,
            const word& phaseName
        )
        :
            basicThermo::implementation(dict, mesh, phaseName),
            rhoThermo::implementation(dict, mesh, phaseName),
            fluidThermo::implementation(dict, mesh, phaseName)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
