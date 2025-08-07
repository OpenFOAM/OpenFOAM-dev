/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::fv::homogeneousNucleation

Description
    Base class for homogeneous nucleation models based on classical nucleation
    theory

SourceFiles
    homogeneousNucleation.C

\*---------------------------------------------------------------------------*/

#ifndef homogeneousNucleation_H
#define homogeneousNucleation_H

#include "phaseChange.H"
#include "phaseSystem.H"
#include "nucleation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class homogeneousNucleation Declaration
\*---------------------------------------------------------------------------*/

class homogeneousNucleation
:
    public phaseChange,
    public nucleation
{
private:

    // Private Data

        //- Phase system
        const phaseSystem& fluid_;

        //- Reference to the solution phase
        const phaseModel& solution_;

        //- Reference to the nucleating phase
        const phaseModel& nucleate_;

        //- Reference to the field associated with the pressure equation
        const volScalarField& p_rgh_;

        //- Counter for the evaluations of the pressure equation sources
        mutable label pressureEquationIndex_;

        //- Diameter of nucleated clusters
        mutable volScalarField::Internal d_;

        //- Phase change rate, per unit volume of the solution
        mutable volScalarField::Internal mDotByAlphaSolution_;

        //- Phase change rate
        mutable volScalarField::Internal mDot_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Correct the stored diameter and phase change rates
        void correctDAndMDot() const;

        //- Calculate and return the diameter of nucleates and the phase change
        //  rate per quantity of the solution phase
        virtual Pair<tmp<volScalarField::Internal>>
            dAndMDotByAlphaSolution() const = 0;


protected:

    // Protected Member Functions

        //- Return the surface tension coefficient between the phases
        tmp<volScalarField::Internal> sigma() const;


public:

    //- Runtime type information
    TypeName("homogeneousNucleation");


    // Constructors

        //- Construct from explicit source name and mesh
        homogeneousNucleation
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // Sources

            //- Return the diameter of nuclei
            virtual tmp<DimensionedField<scalar, volMesh>> d() const;

            //- Return the number rate at which nuclei are generated
            virtual tmp<DimensionedField<scalar, volMesh>> nDot() const;

            //- Return the mass transfer rate
            virtual tmp<DimensionedField<scalar, volMesh>> mDot() const;

            //- Return the nucleation time scale
            virtual tmp<DimensionedField<scalar, volMesh>> tau() const;

            //- Use phaseChange's source functions
            using phaseChange::addSup;

            //- Override the compressible continuity equation to add
            //  linearisation w.r.t alpha
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn
            ) const;


        //- Correct the fvModel
        //  e.g. solve equations, update model, for film, Lagrangian etc.
        virtual void correct();


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
