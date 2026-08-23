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
    Foam::fluidMulticomponentSurfaceThermo

Description
    Base-class for fluid multicomponent Lagrangian surface thermodynamic models

SourceFiles
    fluidMulticomponentSurfaceThermo.C
    fluidMulticomponentSurfaceThermos.C

\*---------------------------------------------------------------------------*/

#ifndef fluidMulticomponentSurfaceThermo_H
#define fluidMulticomponentSurfaceThermo_H

#include "fluidSurfaceThermo.H"
#include "multicomponentSurfaceThermo.H"
#include "SurfaceThermo.H"
#include "MulticomponentSurfaceThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class fluidMulticomponentSurfaceThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidMulticomponentSurfaceThermo
:
    virtual public fluidSurfaceThermo,
    virtual public multicomponentSurfaceThermo
{
protected:

    // Protected Member Functions

        //- Reset to the far thermo
        virtual void reset() const;


public:

    // Public Classes

        //- Forward declare the composite class
        class composite;


    // Public Typedefs

        //- The derived type
        template<class MixtureType>
        using DerivedThermoType =
            MulticomponentSurfaceThermo
            <
                FluidSurfaceThermo
                <
                    SurfaceThermo
                    <
                        MixtureType,
                        composite
                    >
                >
            >;


    //- Runtime type information
    TypeName("fluidMulticomponentSurfaceThermo");


    //- Destructor
    virtual ~fluidMulticomponentSurfaceThermo();
};


/*---------------------------------------------------------------------------*\
            Class fluidMulticomponentSurfaceThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class fluidMulticomponentSurfaceThermo::composite
:
    public surfaceThermo::implementation,
    public fluidSurfaceThermo::implementation,
    public multicomponentSurfaceThermo::implementation,
    public fluidMulticomponentSurfaceThermo
{
public:

    // Constructors

        //- Construct from a reference to the far thermo
        composite(const basicLagrangianSubThermo& farThermo)
        :
            surfaceThermo::implementation(farThermo),
            fluidSurfaceThermo::implementation(farThermo),
            multicomponentSurfaceThermo::implementation(farThermo)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
