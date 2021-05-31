/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::fv::clouds

Description
    Lagrangian clouds fvModel

Usage
    Example usage:
    \verbatim
    clouds
    {
        type    clouds;
    }
    \endverbatim

SourceFiles
    clouds.C

\*---------------------------------------------------------------------------*/

#ifndef clouds_H
#define clouds_H

#include "fvModel.H"
#include "fluidThermo.H"
#include "uniformDimensionedFields.H"
#include "parcelCloudList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class clouds Declaration
\*---------------------------------------------------------------------------*/

class clouds
:
    public fvModel
{
    // Private Data

        //- Reference to the carrier phase thermo
        const fluidThermo& carrierThermo_;

        //- The Lagrangian cloud list
        mutable parcelCloudList clouds_;

        //- Current time index (used for updating)
        mutable label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("clouds");


    // Constructors

        //- Construct from explicit source name and mesh
        clouds
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        clouds
        (
            const clouds&
        ) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Correct

            //- Solve the Lagrangian clouds and update the sources
            virtual void correct();

            //- Prepare for mesh update
            virtual void preUpdateMesh();


        // Add explicit and implicit contributions to compressible equation

            //- Add source to continuity equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

            //- Add source to pressure or enthalpy equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

            //- Add source to momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const clouds&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
