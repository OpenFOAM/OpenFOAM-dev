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

#include "SolidThermo.H"
#include "pureThermo.H"
#include "rhoThermo.H"
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
    virtual public basicThermo,
    virtual public rhoThermo
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
            SolidThermo<BasicThermo<MixtureType, composite>>;

        //- The derived name
        static const word derivedThermoName;


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


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<solidThermo> New
        (
            const fvMesh&,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~solidThermo();


    // Member Functions

        //- Should the dpdt term be included in the enthalpy equation
        virtual Switch dpdt() const = 0;

        //- Return true if thermal conductivity is isotropic
        virtual bool isotropic() const = 0;

        //- Anisotropic thermal conductivity [W/m/K]
        virtual const volVectorField& Kappa() const = 0;
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


public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        implementation(const dictionary&, const fvMesh&, const word&);


    //- Destructor
    virtual ~implementation();


    // Member Functions

        //- The dpdt term should not be included in the enthalpy equation
        virtual Switch dpdt() const
        {
            return false;
        }
};


/*---------------------------------------------------------------------------*\
                   Class solidThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class solidThermo::composite
:
    public basicThermo::implementation,
    public pureThermo,
    public rhoThermo::implementation,
    public solidThermo::implementation
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
            solidThermo::implementation(dict, mesh, phaseName)
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
            composite(dict, mesh, phaseName)
        {}


    // Member Functions

        //- The dpdt term should not be included in the enthalpy equation
        virtual Switch dpdt() const
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
