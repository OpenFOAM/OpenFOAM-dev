/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::diameterModel

Description
    Abstract base-class for dispersed-phase particle diameter models.

SourceFiles
    diameterModel.C
    diameterModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef diameterModel_H
#define diameterModel_H

#include "dictionary.H"
#include "phaseModel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class diameterModel Declaration
\*---------------------------------------------------------------------------*/

class diameterModel
{
    // Private Data

        //- The phase diameter properties dictionary
        dictionary diameterProperties_;

        //- The phase that this model applies
        const phaseModel& phase_;


public:

    //- Runtime type information
    TypeName("diameterModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            diameterModel,
            dictionary,
            (
                const dictionary& diameterProperties,
                const phaseModel& phase
            ),
            (diameterProperties, phase)
        );


    // Constructors

        //- Construct from dictionary and phase
        diameterModel
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~diameterModel();


    // Selectors

        //- Select from dictionary and phase
        static autoPtr<diameterModel> New
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    // Member Functions

        //- Return the phase diameter properties dictionary
        const dictionary& diameterProperties() const
        {
            return diameterProperties_;
        }

        //- Return the phase
        const phaseModel& phase() const
        {
            return phase_;
        }

        //- Return the diameter
        virtual tmp<volScalarField> d() const = 0;

        //- Return the surface area per unit volume
        virtual tmp<volScalarField> a() const = 0;

        //- Correct the model
        virtual void correct();

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
