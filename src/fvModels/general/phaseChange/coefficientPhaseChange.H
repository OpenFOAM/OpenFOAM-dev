/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2025 OpenFOAM Foundation
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
    Foam::fv::coefficientPhaseChange

Description
    This simple model generates a phase change between two phases calculated
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
    coefficientPhaseChange
    {
        type            coefficientPhaseChange;

        phases          (liquid vapour);

        C               [kg/m^2/s] 0.1;
    }
    \endverbatim

SourceFiles
    coefficientPhaseChange.C

\*---------------------------------------------------------------------------*/

#ifndef coefficientPhaseChange_H
#define coefficientPhaseChange_H

#include "phaseChange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class coefficientPhaseChange Declaration
\*---------------------------------------------------------------------------*/

class coefficientPhaseChange
:
    public phaseChange
{
private:

    // Private Data

        //- Phase change coefficient. Multiplies the estimated interfacial
        //  area density return the phase change rate. Units of [kg/m^2/s].
        dimensionedScalar C_;

        //- The source phase fraction field
        const volScalarField& alpha1_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Multiply the argument by the mass fraction of transferring species
        //  in the source phase
        tmp<DimensionedField<scalar, volMesh>> timesY1
        (
            tmp<DimensionedField<scalar, volMesh>> mDot
        ) const;

        //- Return the phase change rate with the source phase and mass
        //  fraction factored out
        virtual tmp<DimensionedField<scalar, volMesh>> mDotByAlpha1Y1() const;


public:

    //- Runtime type information
    TypeName("coefficientPhaseChange");


    // Constructors

        //- Construct from explicit source name and mesh
        coefficientPhaseChange
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // Sources

            //- Return the total phase change rate
            virtual tmp<DimensionedField<scalar, volMesh>> mDot() const;

            //- Return the phase change rate of a specie
            virtual tmp<DimensionedField<scalar, volMesh>> mDot
            (
                const label mDoti
            ) const;

            //- Override the compressible continuity equation to add
            //  linearisation w.r.t alpha
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn
            ) const;

            //- Override the compressible mass fraction equation to add
            //  linearisation w.r.t the mass fraction
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& Yi,
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
