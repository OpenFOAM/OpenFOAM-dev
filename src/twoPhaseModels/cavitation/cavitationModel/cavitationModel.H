/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::cavitationModel

Description
    Abstract base class for cavitation models

SourceFiles
    cavitationModel.C

\*---------------------------------------------------------------------------*/

#ifndef cavitationModel_H
#define cavitationModel_H

#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "fvMatricesFwd.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                                 Class cavitationModel
\*---------------------------------------------------------------------------*/

class cavitationModel
{
protected:

    // Protected data

        //- Mixture properties
        const immiscibleIncompressibleTwoPhaseMixture& mixture_;

        //- Saturation vapour pressure
        dimensionedScalar pSat_;


public:

    //- Runtime type information
    TypeName("cavitationModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            cavitationModel,
            dictionary,
            (
                const dictionary& dict,
                const immiscibleIncompressibleTwoPhaseMixture& mixture
            ),
            (dict, mixture)
        );


    // Constructors

        //- Construct for mixture
        cavitationModel
        (
            const dictionary& dict,
            const immiscibleIncompressibleTwoPhaseMixture& mixture
        );


    // Selector
    static autoPtr<cavitationModel> New
    (
        const dictionary& dict,
        const immiscibleIncompressibleTwoPhaseMixture& mixture
    );


    //- Destructor
    virtual ~cavitationModel()
    {}


    // Member Functions

        //- Return const-access to the saturation vapour pressure
        const dimensionedScalar& pSat() const
        {
            return pSat_;
        }

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField::Internal>> mDotAlphal() const = 0;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (p - pSat)
        virtual Pair<tmp<volScalarField::Internal>> mDotP() const = 0;

        //- Correct the cavitation model
        virtual void correct() = 0;

        //- Read the phaseProperties dictionary and update
        virtual bool read(const dictionary& dict) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
