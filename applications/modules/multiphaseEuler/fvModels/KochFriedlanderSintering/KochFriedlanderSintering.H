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
    Foam::fv::KochFriedlanderSintering

Description
    Sintering model of Koch and Friedlander (1990). The characteristic time for
    sintering is given by

    \f[
        \tau = c_s d_{p_i}^n T^m \exp(T_a/T \cdot [1 - d_{p,min}/d_{p_i}])\;.
    \f]

    Note that the correction factor in the exponential function can be
    eliminated by setting \f$d_{p,min}\f$ to zero which is done by default.

    Reference:
    \verbatim
        Koch, W., & Friedlander, S. K. (1990).
        The effect of particle coalescence on the surface area of a coagulating
        aerosol.
        Journal of Colloid and Interface Science, 140(2), 419-427.
    \endverbatim

Usage
    \table
        Property | Description                | Required    | Default value
        Cs       | Sintering time coefficient | yes         | none
        n        | Particle diameter exponent | yes         | none
        m        | Temperature exponent       | yes         | none
        Ta       | Activation temperature     | yes         | none
        dpMin    | Minimum primary particle diameter | no   | 0
    \endtable

    Example usage:
    \verbatim
    sintering
    {
        type            KochFriedlanderSintering;
        libs            ("libmultiphaseEulerFvModels.so");

        populationBalance aggregates;

        Cs              8.3e24;
        n               4.0;
        m               1.0;
        Ta              3700.0;
    }
    \endverbatim

SourceFiles
    KochFriedlanderSintering.C

\*---------------------------------------------------------------------------*/

#ifndef KochFriedlanderSintering_H
#define KochFriedlanderSintering_H

#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class KochFriedlanderSintering Declaration
\*---------------------------------------------------------------------------*/

class KochFriedlanderSintering
:
    public fvModel
{
private:

    // Private Data

        //- Population balance model
        const diameterModels::populationBalanceModel& popBal_;

        //- Sintering time primary particle diameter exponent
        scalar n_;

        //- Sintering time temperature exponent
        scalar m_;

        //- Sintering time coefficient
        dimensionedScalar Cs_;

        //- Activation temperature
        dimensionedScalar Ta_;

        //- Minimum primary particle diameter, defaults to Zero
        dimensionedScalar dpMin_;

        //- Map from the name of the surface-area-volume-ratio fields to the
        //  index of the size group
        HashTable<label> kappaNameToSizeGroupIndices_;


    // Private Member Functions

        //- Return the dimensions of the sintering time coefficient
        dimensionSet CsDims() const;

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("KochFriedlanderSintering");


    // Constructors

        //- Construct from explicit source name and mesh
        KochFriedlanderSintering
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // Checks

            //- Return true if the fvModel adds a source term to the given
            //  field's transport equation
            virtual bool addsSupToField(const word& fieldName) const;


        // Sources

            //- Return the characteristic time for sintering
            virtual tmp<volScalarField::Internal> tau
            (
                const volScalarField::Internal& kappa
            ) const;

            //- Add a source term to the surface-area-volume-ratio equation
            void addSup
            (
                const volScalarField& alphaFi,
                const volScalarField& rho,
                const volScalarField& kappa,
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
