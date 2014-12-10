/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    A2stract base-class for dispersed-phase particle diameter models.

SourceFiles
    diameterModel.C
    newDiameterModel.C

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

protected:

    // Protected data

        dictionary diameterProperties_;
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

        diameterModel
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~diameterModel();


    // Selectors

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

        //- Return the phase mean diameter field
        virtual tmp<volScalarField> d() const = 0;

        //- Correct the diameter field
        virtual void correct();

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
