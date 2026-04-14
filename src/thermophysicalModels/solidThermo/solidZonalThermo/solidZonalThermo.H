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
    Foam::solidZonalThermo

Description
    Base-class for zonal solid thermodynamic properties.

See also
    Foam::basicThermo

SourceFiles
    solidZonalThermo.C

\*---------------------------------------------------------------------------*/

#ifndef solidZonalThermo_H
#define solidZonalThermo_H

#include "solidThermo.H"
#include "zonalThermo.H"
#include "SolidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class solidZonalThermo Declaration
\*---------------------------------------------------------------------------*/

class solidZonalThermo
:
    virtual public solidThermo,
    virtual public zonalThermo
{
public:

    // Public Classes

        //- Forward declare the composite class
        class composite;


    //- Runtime type information
    TypeName("solidZonalThermo");


    // Public Typedefs

        //- The derived type
        template<class MixtureType>
        using DerivedThermoType =
            SolidThermo<BasicThermo<MixtureType, composite>>;

        //- The derived name
        static word derivedThermoName()
        {
            return "heSolidThermo";
        }


    //- Destructor
    virtual ~solidZonalThermo();
};


/*---------------------------------------------------------------------------*\
                 Class solidZonalThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class solidZonalThermo::composite
:
    public basicThermo::implementation,
    public rhoThermo::implementation,
    public solidThermo::implementation,
    public zonalThermo::implementation,
    public solidZonalThermo
{
public:

    // Constructors

        //- Construct from dictionary, zone names, mesh and phase name
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
            solidThermo::implementation(dict, mesh, phaseName),
            zonalThermo::implementation(dict, mixture.zones(), mesh, phaseName)
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
