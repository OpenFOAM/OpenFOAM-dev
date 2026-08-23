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
    Foam::solidSurfaceThermo

Description
    Base-class for solid Lagrangian surface thermodynamic models

SourceFiles
    solidSurfaceThermo.C
    solidSurfaceThermos.C

\*---------------------------------------------------------------------------*/

#ifndef solidSurfaceThermo_H
#define solidSurfaceThermo_H

#include "surfaceThermo.H"
#include "pureSurfaceThermo.H"
#include "uniformDimensionedFields.H"
#include "uniformGeometricFields.H"
#include "SurfaceThermo.H"
#include "SolidSurfaceThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class solidSurfaceThermo Declaration
\*---------------------------------------------------------------------------*/

class solidSurfaceThermo
:
    virtual public surfaceThermo
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
            SolidSurfaceThermo
            <
                SurfaceThermo
                <
                    MixtureType,
                    composite
                >
            >;


    //- Runtime type information
    TypeName("solidSurfaceThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidSurfaceThermo,
        cloud,
        (const basicLagrangianSubThermo& farThermo),
        (farThermo)
    );


    // Selectors

        //- Select a surface thermo
        static autoPtr<solidSurfaceThermo> New
        (
            const basicLagrangianSubThermo& farThermo
        );


    //- Destructor
    virtual ~solidSurfaceThermo();
};

/*---------------------------------------------------------------------------*\
             Class solidSurfaceThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class solidSurfaceThermo::implementation
:
    virtual public solidSurfaceThermo
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

        //- Construct from a reference to the far thermo
        implementation(const basicLagrangianSubThermo& farThermo);

        //- Disallow default bitwise copy construction
        implementation(const implementation&) = delete;


    //- Destructor
    virtual ~implementation();
};


/*---------------------------------------------------------------------------*\
                 Class solidSurfaceThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class solidSurfaceThermo::composite
:
    public surfaceThermo::implementation,
    public solidSurfaceThermo::implementation,
    public pureSurfaceThermo
{
public:

    // Constructors

        //- Construct from a reference to the far thermo
        composite(const basicLagrangianSubThermo& farThermo)
        :
            surfaceThermo::implementation(farThermo),
            solidSurfaceThermo::implementation(farThermo)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
