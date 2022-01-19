/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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
    Foam::chemistryTabulationMethod

Description
    An abstract class for chemistry tabulation.

SourceFiles
    chemistryTabulationMethod.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryTabulationMethod_H
#define chemistryTabulationMethod_H

#include "IOdictionary.H"
#include "scalarField.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class odeChemistryModel;

/*---------------------------------------------------------------------------*\
                  Class chemistryTabulationMethod Declaration
\*---------------------------------------------------------------------------*/

class chemistryTabulationMethod
{
public:

    //- Runtime type information
    TypeName("chemistryTabulationMethod");


    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        chemistryTabulationMethod,
        dictionary,
        (
            const dictionary& dict,
            const odeChemistryModel& chemistry
        ),
        (dict, chemistry)
    );


    // Constructors

        //- Construct from components
        chemistryTabulationMethod
        (
            const dictionary& dict,
            const odeChemistryModel& chemistry
        );


    // Selectors

        static autoPtr<chemistryTabulationMethod> New
        (
            const IOdictionary& dict,
            const odeChemistryModel& chemistry
        );


    //- Destructor
    virtual ~chemistryTabulationMethod();


    // Member Functions

        //- Return true if the tabulation method tabulates
        //  otherwise return false
        virtual bool tabulates() = 0;

        // Retrieve function: (only virtual here)
        // Try to retrieve a stored point close enough (according to tolerance)
        // to a stored point. If successful, it returns true and store the
        // results in RphiQ, i.e. the result of the integration of phiQ
        virtual bool retrieve
        (
             const scalarField& phiQ,
             scalarField& RphiQ
        ) = 0;

        // Add function: (only virtual here)
        // Add information to the tabulation algorithm. Give the reference for
        // future retrieve (phiQ) and the corresponding result (RphiQ).
        virtual label add
        (
            const scalarField& phiQ,
            const scalarField& RphiQ,
            const label nActive,
            const label li,
            const scalar deltaT
        ) = 0;

        // Update function: (only virtual here)
        // The underlying structure of the tabulation is updated/cleaned
        // to increase the performance of the retrieve
        virtual bool update() = 0;

        virtual void reset() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
