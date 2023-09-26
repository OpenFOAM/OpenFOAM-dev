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
    Foam::fv::filmVoFTransfer

Description
    Film<->VoF transfer model

Usage
    Example usage:
    \verbatim
    filmVoFTransfer
    {
        type    filmVoFTransfer;

        deltaFactorToVoF  1.5;
        alphaToVoF        0.9;

        transferRateCoeff 0.1;
    }
    \endverbatim

SourceFiles
    filmVoFTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef filmVoFTransfer_H
#define filmVoFTransfer_H

#include "fvModel.H"
#include "isothermalFilm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

class VoFFilmTransfer;

/*---------------------------------------------------------------------------*\
                       Class filmVoFTransfer Declaration
\*---------------------------------------------------------------------------*/

class filmVoFTransfer
:
    public fvModel
{
    // Private Data

        //- The film model
        const solvers::isothermalFilm& film_;

        //- Current time index (used for updating)
        mutable label curTimeIndex_;

        //- Factor of the cell height above which the film is transferred
        //  to the VoF
        scalar deltaFactorToVoF_;

        //- VoF limit above which all of the film is transferred to the VoF
        scalar alphaToVoF_;

        //- Transfer rate coefficient
        scalar transferRateCoeff_;

        //- Cached transfer rate
        volScalarField::Internal transferRate_;


    // Private Member Functions

        const VoFFilmTransfer& VoFFilm(const fvModels&) const;

        //- Return the transfer rate from the VoF transferRateFunc
        template<class Type, class TransferRateFunc>
        inline tmp<VolInternalField<Type>> VoFToFilmTransferRate
        (
            TransferRateFunc transferRateFunc,
            const dimensionSet& dimProp
        ) const;

        //- Return the transfer rate of field f
        template<class Type, class FieldType>
        inline tmp<Field<Type>> TransferRate(const FieldType& f) const;


public:

    //- Runtime type information
    TypeName("filmVoFTransfer");


    // Constructors

        //- Construct from explicit source name and mesh
        filmVoFTransfer
        (
            const word& sourceName,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        filmVoFTransfer
        (
            const filmVoFTransfer&
        ) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Correct

            //- Solve the film and update the sources
            virtual void correct();


        // Sources

            //- Add source to phase continuity equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& alpha,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add source to phase energy equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add source to mixture momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volVectorField& U,
                fvMatrix<vector>& eqn
            ) const;


        // Transfer to VoF

            //- Return the volume transfer rate
            tmp<scalarField> transferRate() const;

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
        void operator=(const filmVoFTransfer&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
