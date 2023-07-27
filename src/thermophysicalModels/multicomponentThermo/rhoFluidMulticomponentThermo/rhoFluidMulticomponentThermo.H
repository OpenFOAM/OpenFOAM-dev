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
    Foam::rhoFluidMulticomponentThermo

Description
    Base-class for multi-component fluid thermodynamic properties based on
    density.

See also
    Foam::basicThermo

SourceFiles
    rhoFluidMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef rhoFluidMulticomponentThermo_H
#define rhoFluidMulticomponentThermo_H

#include "rhoFluidThermo.H"
#include "fluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class rhoFluidMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

class rhoFluidMulticomponentThermo
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
        using DerivedThermoType =
            FluidMulticomponentThermo
            <
                MulticomponentThermo
                <
                    RhoFluidThermo
                    <
                        BasicThermo<MixtureType, composite>
                    >
                >
            >;


    //- Runtime type information
    TypeName("rhoFluidMulticomponentThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoFluidMulticomponentThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<rhoFluidMulticomponentThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~rhoFluidMulticomponentThermo();
};


/*---------------------------------------------------------------------------*\
                Class rhoFluidMulticomponentThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class rhoFluidMulticomponentThermo::composite
:
    public basicThermo::implementation,
    public rhoThermo::implementation,
    public fluidThermo::implementation,
    public multicomponentThermo::implementation,
    public rhoFluidMulticomponentThermo
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
                phaseName
            )
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
