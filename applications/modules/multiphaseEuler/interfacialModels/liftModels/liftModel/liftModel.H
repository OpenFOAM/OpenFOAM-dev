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
    Foam::liftModel

Description
    Model for the lift force between two phases

SourceFiles
    liftModel.C
    liftModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef liftModel_H
#define liftModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "BlendedInterfacialModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class liftModel Declaration
\*---------------------------------------------------------------------------*/

class liftModel
{
public:

    //- Runtime type information
    TypeName("liftModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            liftModel,
            dictionary,
            (
                const dictionary& dict,
                const phaseInterface& interface
            ),
            (dict, interface)
        );


    // Static Data Members

        //- Force dimensions
        static const dimensionSet dimF;

        //- Does this model require correcting on fixed flux boundaries?
        static const bool correctFixedFluxBCs = true;


    // Constructors

        //- Construct from a dictionary and an interface
        liftModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~liftModel();


    // Selectors

        static autoPtr<liftModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool outer=true
        );


    // Member Functions

        //- Return lift force
        virtual tmp<volVectorField> F() const = 0;

        //- Return face lift force
        virtual tmp<surfaceScalarField> Ff() const = 0;
};


/*---------------------------------------------------------------------------*\
                      Class blendedLiftModel Declaration
\*---------------------------------------------------------------------------*/

class blendedLiftModel
:
    public BlendedInterfacialModel<liftModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using BlendedInterfacialModel<liftModel>::BlendedInterfacialModel;


    // Selectors

        static autoPtr<blendedLiftModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Member Functions

        //- Return lift force
        tmp<volVectorField> F() const;

        //- Return face lift force
        tmp<surfaceScalarField> Ff() const;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
