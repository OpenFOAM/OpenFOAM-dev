/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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
    Foam::kernelShape

Description
    Base class for ignition kernel flame wrinkling \c Xi correction

Usage
    Example usage:
    \verbatim
    kernelShape
    {
        type            <type>;
    }
    \endverbatim

SourceFiles
    kernelShape.C

\*---------------------------------------------------------------------------*/

#ifndef kernelShape_H
#define kernelShape_H

#include "dimensionedTypes.H"
#include "fvMesh.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class kernelShape Declaration
\*---------------------------------------------------------------------------*/

class kernelShape
{

protected:

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("kernelShape");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            kernelShape,
            dictionary,
            (const fvMesh& mesh, const dictionary& dict),
            (mesh, dict)
        );


    // Constructors

        //- Construct from mesh and dictionary
        kernelShape(const fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        kernelShape(const kernelShape&) = delete;


    // Selectors

        //- Return a reference to the selected kernelShape model
        static autoPtr<kernelShape> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~kernelShape();


    // Member Functions

        //- Return the area of the ignition kernel
        //  calculated from the given kernel volume
        virtual dimensionedScalar Ak(const dimensionedScalar& Vk) const = 0;


    // IO

        //- Update properties from the given XiProperties dictionary
        bool read(const dictionary& XiProperties);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const kernelShape&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
