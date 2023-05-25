/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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
    Foam::relativeVelocityModel

Description

SourceFiles
    relativeVelocityModel.C

\*---------------------------------------------------------------------------*/

#ifndef relativeVelocityModel_H
#define relativeVelocityModel_H

#include "incompressibleDriftFluxMixture.H"
#include "uniformDimensionedFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class relativeVelocityModel Declaration
\*---------------------------------------------------------------------------*/

class relativeVelocityModel
{
    // Private Member Functions

        //- Return the list of patchFieldTypes for Udm derived from U
        wordList UdmPatchFieldTypes() const;


protected:

    // Protected data

        //- Mixture properties
        const incompressibleDriftFluxMixture& mixture_;

        //- Acceleration due to gravity
        const uniformDimensionedVectorField& g_;

        //- Dispersed diffusion velocity
        mutable volVectorField Udm_;


public:

    //- Runtime type information
    TypeName("relativeVelocityModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        relativeVelocityModel,
        dictionary,
        (
            const dictionary& dict,
            const incompressibleDriftFluxMixture& mixture,
            const uniformDimensionedVectorField& g
        ),
        (dict, mixture, g)
    );


    // Constructors

        //- Construct from components
        relativeVelocityModel
        (
            const dictionary& dict,
            const incompressibleDriftFluxMixture& mixture,
            const uniformDimensionedVectorField& g
        );

        //- Disallow default bitwise copy construction
        relativeVelocityModel(const relativeVelocityModel&) = delete;


    // Selector
    static autoPtr<relativeVelocityModel> New
    (
        const dictionary& dict,
        const incompressibleDriftFluxMixture& mixture,
        const uniformDimensionedVectorField& g
    );


    //- Destructor
    virtual ~relativeVelocityModel();


    // Member Functions

        //- Mixture properties
        const incompressibleDriftFluxMixture& mixture() const
        {
            return mixture_;
        }

        //- Return the diffusion velocity of the dispersed phase
        const volVectorField& Udm() const
        {
            return Udm_;
        }

        //- Return the dispersed phase acceleration
        //  Including buoyancy and flow curvature
        tmp<volVectorField> acceleration() const;

        //- Return the stress tensor due to the phase transport
        tmp<volSymmTensorField> tauDm() const;

        //- Return the div stress tensor due to the phase transport
        tmp<volVectorField> divDevTau() const;

        //- Update the diffusion velocity
        virtual void correct() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const relativeVelocityModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
