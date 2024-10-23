/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::XiProfile

Description
    Base class for flame wrinkling profiles

SourceFiles
    XiProfile.C

\*---------------------------------------------------------------------------*/

#ifndef XiProfile_H
#define XiProfile_H

#include "volFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class XiProfile Declaration
\*---------------------------------------------------------------------------*/

class XiProfile
{

protected:

    // Protected data

        const volScalarField& b_;


    // Protected member functions

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict) = 0;


public:

    //- Runtime type information
    TypeName("XiProfile");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            XiProfile,
            dictionary,
            (
                const dictionary& dict,
                const volScalarField& b
            ),
            (
                dict,
                b
            )
        );


    // Constructors

        //- Construct from components
        XiProfile(const volScalarField& b);

        //- Disallow default bitwise copy construction
        XiProfile(const XiProfile&) = delete;


    // Selectors

        //- Return a reference to the selected XiEq model
        static autoPtr<XiProfile> New
        (
            const dictionary& XiProperties,
            const volScalarField& b
        );


    //- Destructor
    virtual ~XiProfile();


    // Member Functions

        //- Return the flame-wrinkling profile
        virtual tmp<volScalarField> profile() const = 0;

        //- Update properties from the given XiProperties dictionary
        bool read(const dictionary& XiProperties);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const XiProfile&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
