/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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
    Foam::fvSource

Description
    Base class for finite volume sources

SourceFiles
    fvSource.C

\*---------------------------------------------------------------------------*/

#ifndef fvSource_H
#define fvSource_H

#include "fvModel.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fvSource Declaration
\*---------------------------------------------------------------------------*/

class fvSource
:
    public fvModel
{
public:

    //- Runtime type information
    TypeName("fvSource");


    // Constructors

        //- Construct from explicit source name and mesh
        fvSource
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        fvSource(const fvSource&) = delete;


    //- Destructor
    virtual ~fvSource();


    // Member Functions

        // Checks

            //- Return true if the fvModel adds a source term to the given
            //  field's transport equation. Must be provided by sources.
            //  Sources potentially apply to every field (or every field of a
            //  phase) and it is not known in advance what all the fields are.
            //  So, this function returns a logical test (e.g., is this the
            //  phase of which there is a source?) rather than testing whether
            //  the name is in a pre-determined set of names.
            virtual bool addsSupToField(const word& fieldName) const = 0;

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation. Sources return null, as all this
            //  logic is handled in addsSupToField.
            virtual wordList addSupFields() const;


        // Access

            //- Return the cells that the source applies to
            virtual labelUList cells() const = 0;

            //- Return the number of cells that the source applies to
            virtual label nCells() const;


        // Sources

            //- Return the source value
            virtual tmp<scalarField> source(const word& fieldName) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
