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
    Foam::fluidMulticomponentCarrierThermo

Description
    Base-class for multicomponent thermodynamic models of the carrier phase of
    a Lagrangian cloud

SourceFiles
    fluidMulticomponentCarrierThermo.C

\*---------------------------------------------------------------------------*/

#ifndef fluidMulticomponentCarrierThermo_H
#define fluidMulticomponentCarrierThermo_H

#include "fluidCarrierThermo.H"
#include "multicomponentCarrierThermo.H"
#include "CarrierThermo.H"
#include "MulticomponentCarrierThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class fluidMulticomponentCarrierThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidMulticomponentCarrierThermo
:
    virtual public fluidCarrierThermo,
    virtual public multicomponentCarrierThermo
{
public:

    // Public Classes

        //- Forward declare the composite class
        class composite;


    // Public Typedefs

        //- The derived type
        template<class MixtureType>
        using DerivedThermoType =
            MulticomponentCarrierThermo
            <
                FluidCarrierThermo
                <
                    CarrierThermo
                    <
                        MixtureType,
                        composite
                    >
                >
            >;


    //- Runtime type information
    TypeName("fluidMulticomponentCarrierThermo");


    //- Destructor
    virtual ~fluidMulticomponentCarrierThermo();
};


/*---------------------------------------------------------------------------*\
            Class fluidMulticomponentCarrierThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class fluidMulticomponentCarrierThermo::composite
:
    public carrierThermo::implementation,
    public fluidCarrierThermo::implementation,
    public multicomponentCarrierThermo::implementation,
    public fluidMulticomponentCarrierThermo
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
            fluidCarrierThermo::implementation(c, carriedCloud, thermo),
            multicomponentCarrierThermo::implementation(c, carriedCloud, thermo)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
