/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2025 OpenFOAM Foundation
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
    Foam::diffusiveMassTransferModel

Description
    Model for diffusive mass transfer coefficients between two phases

SourceFiles
    diffusiveMassTransferModel.C

\*---------------------------------------------------------------------------*/

#ifndef diffusiveMassTransferModel_H
#define diffusiveMassTransferModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "BlendedInterfacialModel.H"
#include "SidedInterfacialModel.H"

namespace Foam
{

class phaseSystem;

/*---------------------------------------------------------------------------*\
                 Class diffusiveMassTransferModel Declaration
\*---------------------------------------------------------------------------*/

class diffusiveMassTransferModel
{
public:

    //- Runtime type information
    TypeName("diffusiveMassTransferModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            diffusiveMassTransferModel,
            dictionary,
            (
                const dictionary& dict,
                const phaseInterface& interface
            ),
            (dict, interface)
        );


    // Static Data Members

        //- Coefficient dimensions
        static const dimensionSet dimK;


    // Constructors

        //- Construct from a dictionary and an interface
        diffusiveMassTransferModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~diffusiveMassTransferModel();


    // Selectors

        static autoPtr<diffusiveMassTransferModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Member Functions

        //- The implicit mass transfer coefficient
        //  Note: this has had the species mass diffusivity factored out
        virtual tmp<volScalarField> K() const = 0;
};


/*---------------------------------------------------------------------------*\
               Class blendedDiffusiveMassTransferModel Declaration
\*---------------------------------------------------------------------------*/

class blendedDiffusiveMassTransferModel
:
    public BlendedInterfacialModel<diffusiveMassTransferModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using
            BlendedInterfacialModel<diffusiveMassTransferModel>::
            BlendedInterfacialModel;


    // Selectors

        inline static autoPtr<blendedDiffusiveMassTransferModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const dictionary& blendingDict
        )
        {
            return autoPtr<blendedDiffusiveMassTransferModel>
            (
                new blendedDiffusiveMassTransferModel
                (
                    dict,
                    interface,
                    blendingDict
                )
            );
        }


    // Member Functions

        //- Return the mass transfer coefficient
        tmp<volScalarField> K() const;
};


/*---------------------------------------------------------------------------*\
            Class sidedBlendedDiffusiveMassTransferModel Declaration
\*---------------------------------------------------------------------------*/

class sidedBlendedDiffusiveMassTransferModel
:
    public SidedInterfacialModel<blendedDiffusiveMassTransferModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using
            SidedInterfacialModel<blendedDiffusiveMassTransferModel>::
            SidedInterfacialModel;


    // Selectors

        inline static autoPtr<sidedBlendedDiffusiveMassTransferModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const dictionary& blendingDict
        )
        {
            return autoPtr<sidedBlendedDiffusiveMassTransferModel>
            (
                new sidedBlendedDiffusiveMassTransferModel
                (
                    dict,
                    interface,
                    blendingDict
                )
            );
        }
};


/*---------------------------------------------------------------------------*\
               Class sidedDiffusiveMassTransferModel Declaration
\*---------------------------------------------------------------------------*/

class sidedDiffusiveMassTransferModel
:
    public SidedInterfacialModel<diffusiveMassTransferModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using
            SidedInterfacialModel<diffusiveMassTransferModel>::
            SidedInterfacialModel;


    // Selectors

        inline static autoPtr<sidedDiffusiveMassTransferModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        )
        {
            return autoPtr<sidedDiffusiveMassTransferModel>
            (
                new sidedDiffusiveMassTransferModel(dict, interface)
            );
        }


    // Member Functions

        //- Return the mass transfer coefficient
        tmp<volScalarField> KinThe(const phaseModel& phase) const;
};


/*---------------------------------------------------------------------------*\
            Class blendedSidedDiffusiveMassTransferModel Declaration
\*---------------------------------------------------------------------------*/

class blendedSidedDiffusiveMassTransferModel
:
    public BlendedInterfacialModel<sidedDiffusiveMassTransferModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using
            BlendedInterfacialModel<sidedDiffusiveMassTransferModel>::
            BlendedInterfacialModel;


    // Selectors

        inline static autoPtr<blendedSidedDiffusiveMassTransferModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const dictionary& blendingDict
        )
        {
            return autoPtr<blendedSidedDiffusiveMassTransferModel>
            (
                new blendedSidedDiffusiveMassTransferModel
                (
                    dict,
                    interface,
                    blendingDict
                )
            );
        }


    // Member Functions

        //- Return the mass transfer coefficient
        tmp<volScalarField> KinThe(const phaseModel& phase) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
