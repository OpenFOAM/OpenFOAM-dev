/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::solidThermo

Description
    Base-class for solid thermodynamic properties.

See also
    Foam::basicThermo

SourceFiles
    solidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef solidThermo_H
#define solidThermo_H

#include "basicThermo.H"
#include "uniformGeometricFields.H"
#include "fvScalarMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class solidThermo Declaration
\*---------------------------------------------------------------------------*/

class solidThermo
:
    virtual public basicThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;

        //- Forward declare the composite class
        class composite;


    //- Runtime type information
    TypeName("solidThermo");


    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidThermo,
        dictionary,
        (const fvMesh& mesh, const dictionary& dict, const word& phaseName),
        (mesh, dict, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<solidThermo> New
        (
            const fvMesh&,
            const word& phaseName = word::null
        );

        //- Selection using a specified dictionary
        static autoPtr<solidThermo> New
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~solidThermo();


    // Member Functions

        //- Density [kg/m^3]
        virtual tmp<volScalarField> rho() const = 0;

        //- Density for patch [kg/m^3]
        virtual tmp<scalarField> rho(const label patchi) const = 0;

        //- Return non-const access to the local density field [kg/m^3]
        virtual volScalarField& rho() = 0;

        //- Old-time density [kg/m^3]
        virtual tmp<volScalarField> rho0() const = 0;

        //- Return true if thermal conductivity is isotropic
        virtual bool isotropic() const = 0;

        //- Anisotropic thermal conductivity [W/m/K]
        virtual tmp<volVectorField> Kappa() const = 0;

        //- Anisotropic thermal conductivity for patch [W/m/K]
        virtual tmp<vectorField> Kappa(const label patchi) const = 0;

        //- Anisotropic thermal conductivity for patch
        //  in the local coordinate system [W/m/K]
        virtual tmp<symmTensorField> KappaLocal(const label patchi) const = 0;

        //- Return the heat flux [W]
        virtual tmp<surfaceScalarField> q() const = 0;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const = 0;
};


/*---------------------------------------------------------------------------*\
                 Class solidThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class solidThermo::implementation
:
    virtual public solidThermo
{
protected:

    // Protected data

        // Fields

            //- Pressure [Pa]
            //  Note: This value should never be used. Solid thermo should only
            //  be instantiated with thermo models that do not depend on
            //  pressure. This uniform field takes a value of NaN, so that if
            //  any thermo models that do depend on pressure are used then the
            //  code will exit.
            uniformGeometricScalarField p_;

            //- Density field [kg/m^3]
            volScalarField rho_;


public:

    // Constructors

        //- Construct from mesh and phase name
        implementation(const fvMesh&, const word& phaseName);

        //- Construct from mesh, dictionary and phase name
        implementation
        (
            const fvMesh&,
            const dictionary& dict,
            const word& phaseName
        );


    //- Destructor
    virtual ~implementation();


    // Member Functions

        //- Density [kg/m^3]
        virtual tmp<volScalarField> rho() const;

        //- Density for patch [kg/m^3]
        virtual tmp<scalarField> rho(const label patchi) const;

        //- Return non-const access to the local density field [kg/m^3]
        virtual volScalarField& rho();

        //- Old-time density [kg/m^3]
        virtual tmp<volScalarField> rho0() const;
};


/*---------------------------------------------------------------------------*\
                   Class solidThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class solidThermo::composite
:
    public basicThermo::implementation,
    public solidThermo::implementation
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
            solidThermo::implementation(mesh, phaseName)
        {}

        //- Construct from mesh, dictionary and phase name
        composite
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& phaseName
        )
        :
            basicThermo::implementation(mesh, dict, phaseName),
            solidThermo::implementation(mesh, dict, phaseName)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
