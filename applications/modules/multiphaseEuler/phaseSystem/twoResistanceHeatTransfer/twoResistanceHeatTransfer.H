/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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
    Foam::fv::twoResistanceHeatTransfer

Description
    Model for heat transfer between two phases. Two heat transfer coefficients
    are used to calculate the heat fluxes that result from the temperature
    differences between the phases and the shared interface.

Usage
    Example usage:
    \verbatim
    heatTransfer
    {
        type            twoResistanceHeatTransfer;
        libs            ("libmultiphaseEulerFvModels.so");

        phases          (air water);

        blending
        {
            type            linear;
            minPartlyContinuousAlpha.air 0;
            minFullyContinuousAlpha.air 1;
            minPartlyContinuousAlpha.water 0;
            minFullyContinuousAlpha.water 1;
        }

        models
        {
            air_dispersedIn_water_inThe_air
            {
                type spherical;
                residualAlpha 1e-4;
            }
            air_dispersedIn_water_inThe_water
            {
                type RanzMarshall;
                residualAlpha 1e-4;
            }
            water_dispersedIn_air_inThe_air
            {
                type RanzMarshall;
                residualAlpha 1e-4;
            }
            water_dispersedIn_air_inThe_water
            {
                type spherical;
                residualAlpha 1e-4;
            }
        }
    }
    \endverbatim

SourceFiles
    twoResistanceHeatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef twoResistanceHeatTransfer_H
#define twoResistanceHeatTransfer_H

#include "phaseSystem.H"
#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
              Class twoResistanceHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class twoResistanceHeatTransfer
:
    public fvModel
{
private:

    // Private Type Definitions

        //- Table of heat transfer models
        typedef
            HashPtrTable
            <
                blendedSidedHeatTransferModel,
                phaseInterfaceKey,
                phaseInterfaceKey::hash
            > modelsTable;

        //- Table of heat transfer coefficients
        typedef
            HashTable
            <
                Pair<tmp<volScalarField>>,
                phaseInterfaceKey,
                phaseInterfaceKey::hash
            > KsTable;


    // Private Data

        //- Phase system
        const phaseSystem& fluid_;

        //- Blended-sided heat transfer coefficient models
        modelsTable models_;

        //- Names of the affected energy fields
        wordList heNames_;

        //- Counter for the evaluations of the energy equation sources
        mutable label energyEquationIndex_;

        //- Cached heat transfer coefficients
        mutable KsTable Ks_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Get the heat transfer coefficients
        const KsTable& Ks() const;


public:

    //- Runtime type information
    TypeName("twoResistanceHeatTransfer");


    // Constructors

        //- Construct from explicit source name and mesh
        twoResistanceHeatTransfer
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // Access

            //- Return the heat transfer coefficients between a pair of phases
            Pair<tmp<volScalarField>> Ks
            (
                const phaseModel& phase1,
                const phaseModel& phase2
            ) const;

            //- Return the heat transfer coefficients between a pair of phases
            //  with a specified residual volume fraction
            Pair<tmp<volScalarField>> Ks
            (
                const phaseModel& phase1,
                const phaseModel& phase2,
                const scalar residualAlpha
            ) const;


        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Sources

            //- Add a source term to the phase-compressible energy equation
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        //- Correct the fvModel
        virtual void correct();


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
