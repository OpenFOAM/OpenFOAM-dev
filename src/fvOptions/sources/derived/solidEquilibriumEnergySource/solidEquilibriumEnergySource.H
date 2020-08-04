/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::fv::solidEquilibriumEnergySource

Description
    This option adds the thermal inertia of a solid phase into the energy
    equation. It assumes that the solid is in thermal equilibrium with the
    surrounding fluid phase.

    The volume fraction of the solid phase is read from constant/alpha.<phase>,
    and the associated thermophysical properties are specified in
    constant/thermophysicalProperties.<phase>.

Usage
    \table
        Property    | Description                            | Req'd? | Default
        phase       | Name of the solid phase                | yes    |
        field       | Name of the energy field to apply the option to \\
                                                             | yes    |
    \endtable

    Example specification:
    \verbatim
    <fvOptionName>
    {
        type        solidEquilibriumEnergySource;
        phase       solid;
        field       e;
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef solidEquilibriumEnergySource_H
#define solidEquilibriumEnergySource_H

#include "fvOption.H"
#include "volFields.H"
#include "solidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class solidEquilibriumEnergySource Declaration
\*---------------------------------------------------------------------------*/

class solidEquilibriumEnergySource
:
    public option
{
    // Private Member Data

        //- The name of the phase
        const word phaseName_;

        //- Get the volume fraction field
        const volScalarField& alpha() const;

        //- Get the thermo
        const solidThermo& thermo() const;


public:

    //- Runtime type information
    TypeName("solidEquilibriumEnergySource");


    // Constructors

        //- Construct from components
        solidEquilibriumEnergySource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        solidEquilibriumEnergySource
        (
            const solidEquilibriumEnergySource&
        ) = delete;


    //- Destructor
    virtual ~solidEquilibriumEnergySource();


    // Member Functions

        // Evaluation

            //- Explicit and implicit sources for compressible equations
            virtual void addSup
            (
                const volScalarField&,
                fvMatrix<scalar>&,
                const label
            ) const;

            //- Explicit and implicit sources for phase equations
            virtual void addSup
            (
                const volScalarField&,
                const volScalarField&,
                fvMatrix<scalar>&,
                const label
            ) const;


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solidEquilibriumEnergySource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
