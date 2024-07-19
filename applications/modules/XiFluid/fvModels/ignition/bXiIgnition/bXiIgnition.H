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
    Foam::fv::bXiIgnition

Description
    Abstract base-class for ignition models
    for the Weller b-Xi combustion models

SourceFiles
    bXiIgnition.C

\*---------------------------------------------------------------------------*/

#ifndef bXiIgnition_H
#define bXiIgnition_H

#include "fvModel.H"
#include "XiFluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class bXiIgnition Declaration
\*---------------------------------------------------------------------------*/

class bXiIgnition
:
    public fvModel
{
    // Private Data

        //- The XiFluid solver
        const solvers::XiFluid& XiFluid_;


public:

    //- Runtime type information
    TypeName("bXiIgnition");


    // Constructors

        //- Construct from explicit source name and mesh
        bXiIgnition
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        bXiIgnition
        (
            const bXiIgnition&
        ) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;

            virtual bool igniting() const = 0;

            virtual bool ignited() const = 0;


        // Add explicit and implicit contributions to compressible equation

            virtual void XiCorr
            (
                volScalarField& Xi,
                const volScalarField& b,
                const volScalarField& mgb
            ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const bXiIgnition&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
