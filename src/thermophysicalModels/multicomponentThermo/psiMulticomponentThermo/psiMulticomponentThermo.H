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
    Foam::psiMulticomponentThermo

Description
    Base-class for multi-component fluid thermodynamic properties based on
    compressibility.

See also
    Foam::basicThermo

SourceFiles
    psiMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef psiMulticomponentThermo_H
#define psiMulticomponentThermo_H

#include "psiThermo.H"
#include "fluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class psiMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

class psiMulticomponentThermo
:
    virtual public psiThermo,
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
                    PsiThermo
                    <
                        BasicThermo<MixtureType, composite>
                    >
                >
            >;


    //- Runtime type information
    TypeName("psiMulticomponentThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        psiMulticomponentThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<psiMulticomponentThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~psiMulticomponentThermo();
};


/*---------------------------------------------------------------------------*\
                Class psiMulticomponentThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class psiMulticomponentThermo::composite
:
    public basicThermo::implementation,
    public fluidThermo::implementation,
    public psiThermo::implementation,
    public multicomponentThermo::implementation,
    public psiMulticomponentThermo
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
            fluidThermo::implementation(dict, mesh, phaseName),
            psiThermo::implementation(dict, mesh, phaseName),
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
