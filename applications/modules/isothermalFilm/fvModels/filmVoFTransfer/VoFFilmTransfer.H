/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::fv::VoFFilmTransfer

Description
    Film<->VoF transfer model

Usage
    Example usage:
    \verbatim
    VoFFilmTransfer
    {
        type        VoFFilmTransfer;

        libs        ("libfilmVoFTransfer.so");

        filmPatch   film;
        phase       liquid;

        deltaFactorToFilm 0.9;
        alphaToFilm       0.86;

        transferRateCoeff 0.1;
    }
    \endverbatim

SourceFiles
    VoFFilmTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef VoFFilmTransfer_H
#define VoFFilmTransfer_H

#include "fvModel.H"
#include "compressibleVoF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class VoFFilmTransfer Declaration
\*---------------------------------------------------------------------------*/

class VoFFilmTransfer
:
    public fvModel
{
    // Private Data

        //- The VoF model
        const solvers::compressibleVoF& VoF_;

        //- Film patch name
        word filmPatchName_;

        const label filmPatchi_;

        //- The name of the transferred phase
        word phaseName_;

        //- Reference to the transferred phase thermo
        const rhoFluidThermo& thermo_;

        //- Reference to the transferred phase volume fraction
        const volScalarField& alpha_;

        //- Current time index (used for updating)
        mutable label curTimeIndex_;

        //- Factor of the cell height below which the VoF may be transferred
        //  to the film
        scalar deltaFactorToFilm_;

        //- VoF limit below which the VoF may be transferred to the film
        scalar alphaToFilm_;

        //- Transfer rate coefficient
        scalar transferRateCoeff_;

        volScalarField::Internal transferRate_;


    // Private Member Functions

        //- Return the transfer rate from the film transferRateFunc
        template<class Type, class TransferRateFunc>
        inline tmp<VolInternalField<Type>> filmVoFTransferRate
        (
            TransferRateFunc transferRateFunc,
            const dimensionSet& dimProp
        ) const;

        //- Return the transfer rate of field f
        template<class Type, class FieldType>
        inline tmp<Field<Type>> TransferRate(const FieldType& f) const;


public:

    //- Runtime type information
    TypeName("VoFFilmTransfer");


    // Constructors

        //- Construct from explicit source name and mesh
        VoFFilmTransfer
        (
            const word& sourceName,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        VoFFilmTransfer
        (
            const VoFFilmTransfer&
        ) = delete;


    // Member Functions

            label filmPatchIndex() const
            {
                return filmPatchi_;
            }

            const volScalarField& alpha() const
            {
                return alpha_;
            }


        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Correct

            //- Solve the film and update the sources
            virtual void correct();


        // Sources

            //- Add implicit contribution to phase continuity equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add implicit contribution to phase energy equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add implicit contribution to mixture momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volVectorField& U,
                fvMatrix<vector>& eqn
            ) const;


        // Transfer to film

            //- Return the mass transfer rate
            tmp<scalarField> rhoTransferRate() const;

            //- Return the energy transfer rate
            tmp<scalarField> heTransferRate() const;

            //- Return the momentum transfer rate
            tmp<vectorField> UTransferRate() const;


        // Mesh changes

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const VoFFilmTransfer&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
