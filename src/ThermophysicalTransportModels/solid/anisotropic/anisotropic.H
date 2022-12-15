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

    If the coordinate transformed kappa does not align exactly with the boundary
    because the patch face orientations do not conform to the coordinate system
    exactly it may be beneficial for convergence and accuracy to enforce
    alignment at the boundary by setting the optional \c boundaryAligned to
    true.

Usage
    Example of the anisotropic thermal conductivity specification in
    thermophysicalProperties with two zone-based coordinate systems in
    addition to the default:

    \verbatim
    model anisotropic;

    // Force aligned handling of kappa irrespective
    // of the calculated patch alignment factors.
    boundaryAligned true;

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
#include "MeshObjects.H"

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
    public solidThermophysicalTransportModel,
    public UpdateableMeshObject<fvMesh>
{
    // Private member data

        coordinateSystem coordinateSystem_;

        //- Optional zone coordinate systems
        PtrDictionary<coordinateSystem> zoneCoordinateSystems_;

        //- Switch to override boundary alignment checks
        //  and force aligned handling of kappa irrespective
        //  of the calculated patch alignment factors
        Switch boundaryAligned_;

        // List of patch faces adjacent to coordinate zones
        mutable labelListListList zonesPatchFaces_;

        //- Alignment of Kappa for each patch
        //  true is closely aligned, zero otherwise
        boolList aligned_;

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

        //- Return the patch heat flux correction [W/m^2]
        //  For patch-aligned thermal conductivity qCorr is null
        virtual tmp<scalarField> qCorr(const label patchi) const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

        //- Correct the anisotropic viscosity
        virtual void predict();


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap& map);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap& map);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap& map);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidThermophysicalTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
