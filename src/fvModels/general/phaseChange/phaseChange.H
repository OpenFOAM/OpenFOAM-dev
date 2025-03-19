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
    Foam::fv::phaseChange

Description
    Base class for phase change models

SourceFiles
    phaseChange.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChange_H
#define phaseChange_H

#include "hashedWordList.H"
#include "massTransfer.H"
#include "ThermoRefPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class basicThermo;
class fluidThermo;
class multicomponentThermo;
class fluidMulticomponentThermo;

namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class phaseChange Declaration
\*---------------------------------------------------------------------------*/

class phaseChange
:
    public massTransfer
{
private:

    // Private Data

        //- The thermo references
        const ThermoRefPair<basicThermo> thermos_;

        //- Names of the energy fields
        const Pair<word> heNames_;

        //- The names of the transferring species
        hashedWordList species_;

        //- The indices of the transferring species in the two phases, or -1 if
        //  the phase is not multicomponent
        List<labelPair> specieis_;

        //- Whether or not to linearise the energy source
        bool energySemiImplicit_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Initialise the indices of the transferring species in the two
        //  phases, or -1 if the phase is not multicomponent
        const List<labelPair> initSpecieis() const;


protected:

    // Protected Member Functions

        //- Read the names of the transferring specie
        wordList readSpecie(const dictionary& dict, const bool required) const;

        //- Read the names of the transferring species
        wordList readSpecies(const dictionary& dict, const bool required) const;

        //- Re-read the names of the transferring specie
        void reReadSpecie(const dictionary& dict) const;

        //- Re-read the names of the transferring species
        void reReadSpecies(const dictionary& dict) const;

        //- Set the names of the transferring species
        void setSpecies
        (
            const word& name,
            const word& modelType,
            const wordList& species
        );

        //- Set the names of the transferring species
        void setSpecies(const wordList& species);

        //- Re-set the names of the transferring species
        void reSetSpecies(const wordList& species);

        //- Access the pressure field
        const volScalarField& p() const;

        //- Add a boundary field to the given internal field
        static tmp<volScalarField> vifToVf
        (
            const tmp<DimensionedField<scalar, volMesh>>& tvif
        );

        //- Remove the boundary field from the given geometric field
        static tmp<DimensionedField<scalar, volMesh>> vfToVif
        (
            const tmp<volScalarField>& tvf
        );


public:

    //- Runtime type information
    TypeName("phaseChange");


    // Constructors

        //- Construct from explicit source name and mesh
        phaseChange
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const wordList& species
        );


    // Member Functions

        // Access

            //- Return the thermo references
            inline const ThermoRefPair<basicThermo>& thermos() const;

            //- Return the fluid thermo references
            const ThermoRefPair<fluidThermo>
                fluidThermos(const bool, const bool) const;

            //- Return the multicomponent thermo references
            const ThermoRefPair<multicomponentThermo>
                multicomponentThermos(const bool, const bool) const;

            //- Return the fluid multicomponent thermo references
            const ThermoRefPair<fluidMulticomponentThermo>
                fluidMulticomponentThermos(const bool, const bool) const;

            //- Return the names of the energy fields
            inline const Pair<word>& heNames() const;

            //- Return the names of the transferring species. Empty if neither
            //  thermo is multicomponent.
            inline const hashedWordList& species() const;

            //- Return the indices of the transferring specie in the two
            //  phases, or -1 if the phase is not multicomponent
            const labelPair& specieis(const label mDoti = -1) const;


        // Evaluation

            //- Return the temperature at which the phases are considered to be
            //  changing. By default this is considered to be the temperature
            //  of the "source" phase (i.e., the phase for which mDot is
            //  negative), but this can be overridden to account for heat
            //  transfer modelling or similar.
            virtual tmp<DimensionedField<scalar, volMesh>> Tchange() const;

            //- Return the fraction of the latent heat that is transferred into
            //  the second phase. By default this is weighted by the phase
            //  thermal conductivities, but this can be overridden to account
            //  for heat transfer modelling or similar.
            virtual tmp<DimensionedField<scalar, volMesh>> Lfraction() const;

            //- Return the latent heat
            tmp<DimensionedField<scalar, volMesh>> L
            (
                const label mDoti = -1
            ) const;

            //- Return the latent heat for a given changing temperature
            tmp<DimensionedField<scalar, volMesh>> L
            (
                const volScalarField::Internal& Tchange,
                const label mDoti = -1
            ) const;

            //- Return the latent heat for a patch and a given changing
            //  temperature
            tmp<scalarField> L
            (
                const label patchi,
                const scalarField& Tchange,
                const label mDoti = -1
            ) const;


        // Sources

            //- Return the total phase change rate
            virtual tmp<DimensionedField<scalar, volMesh>> mDot() const;

            //- Return the mass transfer rate of a specie
            virtual tmp<DimensionedField<scalar, volMesh>> mDot
            (
                const label mDoti
            ) const;

            //- Override the energy equation to add the phase change heat, or
            //  the species equations to add the relevant mass sources
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& heOrYi,
                fvMatrix<scalar>& eqn
            ) const;


        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseChangeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
