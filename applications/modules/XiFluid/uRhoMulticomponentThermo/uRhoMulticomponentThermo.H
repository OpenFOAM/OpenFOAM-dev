/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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
    Foam::uRhoMulticomponentThermo

Description
    Base-class for combustion fluid thermodynamic properties based on
    compressibility.

See also
    Foam::basicThermo

SourceFiles
    uRhoMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef uRhoMulticomponentThermo_H
#define uRhoMulticomponentThermo_H

#include "rhoFluidThermo.H"
#include "fluidMulticomponentThermo.H"
#include "URhoFluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class uRhoMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

class uRhoMulticomponentThermo
:
    virtual public rhoFluidThermo,
    virtual public fluidMulticomponentThermo
{

public:

    // Public Classes

        //- Forward declare the composite class
        class composite;


    // Public Typedefs

        //- The derived type
        template<class MixtureType>
        using DerivedThermoType = URhoFluidMulticomponentThermo
        <
            FluidMulticomponentThermo
            <
                MulticomponentThermo
                <
                    BasicThermo<MixtureType, composite>
                >
            >
        >;

        //- The derived name
        static const word derivedThermoName()
        {
            return "uRhoThermo";
        }


    //- Runtime type information
    TypeName("uRhoMulticomponentThermo");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            uRhoMulticomponentThermo,
            fvMesh,
            (const fvMesh& mesh, const word& phaseName),
            (mesh, phaseName)
        );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<uRhoMulticomponentThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~uRhoMulticomponentThermo();


    // Member Functions

        //- Return the fuel-oxidant equivalence ratio
        virtual tmp<volScalarField> Phi() const = 0;

        //- Return the prompt burnt gas composition
        virtual PtrList<volScalarField::Internal> prompt() const = 0;

        //- Reset the mixture to an unburnt state and update EGR
        virtual void reset
        (
            volScalarField& b,
            volScalarField& c,
            const PtrList<volScalarField>& Yb,
            const volScalarField& heb
        ) = 0;
};


/*---------------------------------------------------------------------------*\
              Class uRhoMulticomponentThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class uRhoMulticomponentThermo::composite
:
    public basicThermo::implementation,
    public rhoThermo::implementation,
    public fluidThermo::implementation,
    public multicomponentThermo::implementation,
    public uRhoMulticomponentThermo
{
public:

    // Constructors

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
            fluidThermo::implementation(dict, mesh, phaseName),
            multicomponentThermo::implementation
            (
                dict,
                mixture.specieNames(),
                mesh,
                phaseName,
                false
            )
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
