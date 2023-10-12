/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::profileModel

Description
    Base class for profile models

SourceFiles
    profileModel.C

\*---------------------------------------------------------------------------*/

#ifndef profileModel_H
#define profileModel_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class profileModel Declaration
\*---------------------------------------------------------------------------*/

class profileModel
{

protected:

    // Protected data

        //- Coefficients dictionary
        const dictionary dict_;

        //- Name of profile model
        const word name_;

        //- File name (optional)
        fileName fName_;


    // Protected Member Functions

        //- Return true if file name is set
        bool readFromFile() const;


public:

    //- Runtime type information
    TypeName("profileModel");


        // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            profileModel,
            dictionary,
            (
                const dictionary& dict,
                const word& modelName
            ),
            (dict, modelName)
        );


    // Selectors

        //- Return a reference to the selected fvModel model
        static autoPtr<profileModel> New(const dictionary& dict);


    //- Constructor
    profileModel(const dictionary& dict, const word& modelName);


    //- Destructor
    virtual ~profileModel();


    // Member Functions

        // Access

            //- Return const access to the source name
            const word& name() const;


        // Evaluation

            //- Return the Cd and Cl for a given angle-of-attack
            virtual void Cdl
            (
                const scalar alpha,
                scalar& Cd,
                scalar& Cl
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
