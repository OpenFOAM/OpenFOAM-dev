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
    Foam::twoPhaseChangeModel

Description

SourceFiles
    twoPhaseChangeModel.C
    phaseChangeModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseChangeModel_H
#define twoPhaseChangeModel_H

#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "fvMatricesFwd.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class twoPhaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseChangeModel
:
    public IOdictionary
{
    // Private Member Functions

        //- Construct the base IO object
        IOobject createIOobject
        (
            const immiscibleIncompressibleTwoPhaseMixture& mixture
        ) const;


protected:

    // Protected data

        //- Reference to the two-phase mixture
        const immiscibleIncompressibleTwoPhaseMixture& mixture_;

        //- Model coefficient dictionary
        dictionary twoPhaseChangeModelCoeffs_;


public:

    //- Runtime type information
    TypeName("phaseChangeModel");

    //- Default phaseChangeProperties dictionary name
    static const word phaseChangePropertiesName;


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            twoPhaseChangeModel,
            dictionary,
            (
                const immiscibleIncompressibleTwoPhaseMixture& mixture
            ),
            (mixture)
        );


    // Constructors

        //- Construct from dictionary
        twoPhaseChangeModel
        (
            const word& type,
            const immiscibleIncompressibleTwoPhaseMixture& mixture
        );

        //- Disallow default bitwise copy construction
        twoPhaseChangeModel(const twoPhaseChangeModel&);


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<twoPhaseChangeModel> New
        (
            const immiscibleIncompressibleTwoPhaseMixture& mixture
        );


    //- Destructor
    virtual ~twoPhaseChangeModel()
    {}


    // Member Functions

        //- Return the phase-change explicit and implicit sources
        //  for the phase-fraction equation
        virtual Pair<tmp<volScalarField::Internal>> Salpha
        (
            volScalarField& alpha
        ) const = 0;

        //- Return the phase-change source matrix
        //  for the p_rgh pressure equation
        virtual tmp<fvScalarMatrix> Sp_rgh
        (
            const volScalarField& rho,
            const volScalarField& gh,
            volScalarField& p_rgh
        ) const = 0;

        //- Return the phase-change source matrix for the momentum equation
        virtual tmp<fvVectorMatrix> SU
        (
            const volScalarField& rho,
            const surfaceScalarField& rhoPhi,
            volVectorField& U
        ) const = 0;

        //- Correct the phaseChange model
        virtual void correct() = 0;

        //- Read the transportProperties dictionary and update
        virtual bool read() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const twoPhaseChangeModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
