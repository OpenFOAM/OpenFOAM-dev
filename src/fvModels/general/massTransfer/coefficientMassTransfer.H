/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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
    Foam::fv::coefficientMassTransfer

Description
    This simple model generates a mass transfer between two phases calculated
    from the following expression:

        \f[
            \dot{m}/V = C \alpha \grad \alpha
        \f]

    \vartable
        \dot{m}/V | mass transfer rate per unit volume
        C         | coefficient
        \alpha    | volume fraction of the source phase
    \endvartable

Usage
    Example usage:
    \verbatim
    coefficientMassTransfer
    {
        type            coefficientMassTransfer;

        phases          (liquid vapour);

        C               [kg/m^2/s] 0.1;
    }
    \endverbatim

SourceFiles
    coefficientMassTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef coefficientMassTransfer_H
#define coefficientMassTransfer_H

#include "massTransfer.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class coefficientMassTransfer Declaration
\*---------------------------------------------------------------------------*/

class coefficientMassTransfer
:
    public massTransfer
{
private:

    // Private Data

        //- Mass flow rate coefficient. Multiplies the estimated interfacial
        //  area density return the mass transfer rate. Units of [kg/m^2/s].
        dimensionedScalar C_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("coefficientMassTransfer");


    // Constructors

        //- Construct from explicit source name and mesh
        coefficientMassTransfer
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // Sources

            //- Return the mass transfer rate
            virtual tmp<DimensionedField<scalar, volMesh>> mDot() const;

            //- Override the incompressible continuity equation to add
            //  linearisation w.r.t alpha
            void addSup
            (
                const volScalarField& alpha,
                fvMatrix<scalar>& eqn
            ) const;

            //- Override the compressible continuity equation to add
            //  linearisation w.r.t alpha
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn
            ) const;


        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
