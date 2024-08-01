/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::fv::constantbXiIgnition

Description
    Simple constant rate ignition model for the Weller b-Xi combustion models

Usage
    Example usage:
    \verbatim
    constantbXiIgnition1
    {
        type            constantbXiIgnition;

        select          cellZone;
        cellZone        ignition;

        start           0;
        duration        0.003;
        strength        2;

        XiCorr
        {
            type            spherical;
            select          all;
        }
    }
    \endverbatim

    Where:
    \table
        Property     | Description                 | Required   | Default value
        select       | Ignition cell set selection | yes        |
        start        | Ignition start time         | yes        |
        duration     | Ignition duration           | yes        |
        strength     | Ignition strength [1/s]     | yes        |
        XiCorr       | Flame-wrinkling correction  | yes        |
    \endtable

See also
    Foam::XiCorrModel
    Foam::XiCorrModels::planar
    Foam::XiCorrModels::cylindrical
    Foam::XiCorrModels::spherical

SourceFiles
    constantbXiIgnition.C

\*---------------------------------------------------------------------------*/

#ifndef constantbXiIgnition_H
#define constantbXiIgnition_H

#include "bXiIgnition.H"
#include "XiCorrModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class constantbXiIgnition Declaration
\*---------------------------------------------------------------------------*/

class constantbXiIgnition
:
    public bXiIgnition
{
    // Private Data

        //- The set of ignition cells
        fvCellSet set_;

        autoPtr<XiCorrModel> XiCorrModel_;

        //- Ignition start time
        dimensionedScalar start_;

        //- Ignition duration
        dimensionedScalar duration_;

        //- Ignition strength
        dimensionedScalar strength_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("constantbXiIgnition");


    // Constructors

        //- Construct from explicit source name and mesh
        constantbXiIgnition
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        constantbXiIgnition
        (
            const constantbXiIgnition&
        ) = delete;


    // Member Functions

        // Checks

            virtual bool igniting() const;

            virtual bool ignited() const;


        // Add explicit and implicit contributions to compressible equation

            //- Add ignition contribution to b equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& b,
                fvMatrix<scalar>& eqn
            ) const;

            virtual void XiCorr
            (
                volScalarField& Xi,
                const volScalarField& b,
                const volScalarField& mgb
            ) const;


        // Mesh motion

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const constantbXiIgnition&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
