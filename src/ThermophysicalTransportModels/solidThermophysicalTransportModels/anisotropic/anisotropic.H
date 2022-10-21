/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::solidThermophysicalTransportModels::anisotropic

Description
    Solid thermophysical transport model for anisotropic thermal conductivity

    The anisotropic thermal conductivity field is evaluated from the solid
    material anisotropic kappa specified in the physicalProperties dictionary
    transformed into the global coordinate system using default
    coordinate system and optionally additional coordinate systems specified
    per-zone in the thermophysicalProperties dictionary.

Usage
    Example of the anisotropic thermal conductivity specification in
    thermophysicalProperties with two zone-based coordinate systems in
    addition to the default:

    \verbatim
    model anisotropic;

    // Default coordinate system
    coordinateSystem
    {
        type        cartesian;
        origin      (0 0 0);
        coordinateRotation
        {
            type        cylindrical;
            e3          (1 0 0);
        }
    }

    // Optional zone coordinate systems
    zones
    {
        coil1
        {
            type        cartesian;
            origin      (0.1 0.2 0.7);
            coordinateRotation
            {
                type        cylindrical;
                e3          (0.5 0.866 0);
            }
        }

        coil2
        {
            type        cartesian;
            origin      (0.4 0.5 1);
            coordinateRotation
            {
                type        cylindrical;
                e3          (0.866 0.5 0);
            }
        }
    }
    \endverbatim

SourceFiles
    anisotropic.C

\*---------------------------------------------------------------------------*/

#ifndef anisotropic_H
#define anisotropic_H

#include "solidThermophysicalTransportModel.H"
#include "coordinateSystem.H"
#include "PtrDictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
                           Class anisotropic Declaration
\*---------------------------------------------------------------------------*/

class anisotropic
:
    public solidThermophysicalTransportModel
{
    // Private member data

        coordinateSystem coordinateSystem_;

        //- Optional zone coordinate systems
        PtrDictionary<coordinateSystem> zoneCoordinateSystems_;

        // List of patch faces adjacent to coordinate zones
        mutable labelListListList zonesPatchFaces_;

        //- Find all the patch faces adjacent to zones
        void setZonesPatchFaces() const;

        //- Thermal conductivity [W/m/K]
        tmp<volSymmTensorField> Kappa() const;

        //- Thermal conductivity for patch [W/m/K]
        tmp<symmTensorField> Kappa(const label patchi) const;


public:

    //- Runtime type information
    TypeName("anisotropic");


    // Constructors

        //- Construct from solid thermophysical properties
        anisotropic(const solidThermo& thermo);


    //- Destructor
    virtual ~anisotropic()
    {}


    // Member Functions

        //- Read thermophysicalTransport dictionary
        virtual bool read();

        //- Thermal conductivity [W/m/K]
        virtual tmp<volScalarField> kappa() const;

        //- Thermal conductivity for patch [W/m/K]
        virtual tmp<scalarField> kappa(const label patchi) const;

        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

        //- Correct the anisotropic viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidThermophysicalTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
