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
    Foam::bRhoMulticomponentThermo

Description
    Base-class for combustion fluid thermodynamic properties based on
    compressibility.

See also
    Foam::basicThermo

SourceFiles
    bRhoMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef bRhoMulticomponentThermo_H
#define bRhoMulticomponentThermo_H

#include "rhoFluidThermo.H"
#include "fluidMulticomponentThermo.H"
#include "BRhoFluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class bRhoMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

class bRhoMulticomponentThermo
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
        using DerivedThermoType = BRhoFluidMulticomponentThermo
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
            return "bRhoThermo";
        }


    //- Runtime type information
    TypeName("bRhoMulticomponentThermo");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            bRhoMulticomponentThermo,
            fvMesh,
            (const fvMesh& mesh, const word& phaseName),
            (mesh, phaseName)
        );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<bRhoMulticomponentThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~bRhoMulticomponentThermo();


    // Member Functions

        //- Standard enthalpy of formation [J/kg]
        virtual tmp<volScalarField> hf() const = 0;
};


/*---------------------------------------------------------------------------*\
              Class bRhoMulticomponentThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class bRhoMulticomponentThermo::composite
:
    public basicThermo::implementation,
    public rhoThermo::implementation,
    public fluidThermo::implementation,
    public multicomponentThermo::implementation,
    public bRhoMulticomponentThermo
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
