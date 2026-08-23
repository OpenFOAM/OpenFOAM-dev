/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::solidCarrierThermo

Description
    Base-class for multicomponent thermodynamic models of the carrier phase of
    a Lagrangian cloud

SourceFiles
    solidCarrierThermo.C

\*---------------------------------------------------------------------------*/

#ifndef solidCarrierThermo_H
#define solidCarrierThermo_H

#include "carrierThermo.H"
#include "pureCarrierThermo.H"
#include "uniformDimensionedFields.H"
#include "uniformGeometricFields.H"
#include "CarrierThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class solidCarrierThermo Declaration
\*---------------------------------------------------------------------------*/

class solidCarrierThermo
:
    virtual public carrierThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;

        //- Forward declare the composite class
        class composite;


    // Public Typedefs

        //- The derived type
        template<class MixtureType>
        using DerivedThermoType =
            CarrierThermo<MixtureType, composite>;


    //- Runtime type information
    TypeName("solidCarrierThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidCarrierThermo,
        cloud,
        (
            const cloud& c,
            const clouds::carried& carriedCloud,
            const basicThermo& thermo
        ),
        (c, carriedCloud, thermo)
    );


    // Selectors

        //- Select a carrier thermo
        static autoPtr<solidCarrierThermo> New
        (
            const cloud& c,
            const clouds::carried& carriedCloud,
            const word& phaseName
        );


    //- Destructor
    virtual ~solidCarrierThermo();
};


/*---------------------------------------------------------------------------*\
            Class solidCarrierThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class solidCarrierThermo::implementation
:
    public virtual solidCarrierThermo
{
protected:

    // Protected Data

        //- Pressure [Pa]
        //  Note: This value should never be used. Solid thermo should only be
        //  instantiated with thermo models that do not depend on pressure.
        //  This uniform field takes a value of NaN, so that if any thermo
        //  models that do depend on pressure are used then the code will exit.
        uniformGeometricScalarField p_;


    // Protected Member Functions

        //- Pressure for a sub-mesh [Pa]
        virtual tmp<uniformDimensionedScalarField> p
        (
            const LagrangianSubMesh&
        ) const;


public:

    // Constructors

        //- Construct from a reference to the cloud and its base classes
        implementation
        (
            const cloud& c,
            const clouds::carried& carriedCloud,
            const basicThermo& thermo
        );


    //- Destructor
    virtual ~implementation();
};


/*---------------------------------------------------------------------------*\
            Class solidCarrierThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class solidCarrierThermo::composite
:
    public carrierThermo::implementation,
    public solidCarrierThermo::implementation,
    public pureCarrierThermo
{
public:

    // Constructors

        //- Construct from a reference to the cloud and its base classes
        composite
        (
            const cloud& c,
            const clouds::carried& carriedCloud,
            const basicThermo& thermo
        )
        :
            carrierThermo::implementation(c, carriedCloud, thermo),
            solidCarrierThermo::implementation(c, carriedCloud, thermo)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
