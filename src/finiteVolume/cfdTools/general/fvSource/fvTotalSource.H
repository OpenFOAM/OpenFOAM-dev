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
    Foam::fvTotalSource

Description
    Base class for sources which are specified as a total value (e.g., volume
    or mass flow rate), rather than a specific value (e.g., mass flow rate per
    unit volume).

SourceFiles
    fvTotalSource.C

\*---------------------------------------------------------------------------*/

#ifndef fvTotalSource_H
#define fvTotalSource_H

#include "fvSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvCellSet;

/*---------------------------------------------------------------------------*\
                        Class fvTotalSource Declaration
\*---------------------------------------------------------------------------*/

class fvTotalSource
:
    public fvSource
{
private:

    // Private Data

        //- Name of the phase
        word phaseName_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


protected:

    // Protected Member Functions

        // Sources

            //- Add a source term to a field-less proxy equation
            void addSource(fvMatrix<scalar>& eqn) const;

            //- Add a source term to an equation
            template<class Type>
            void addSupType
            (
                const VolField<Type>& field,
                fvMatrix<Type>& eqn
            ) const;

            //- Add a source term to a compressible equation
            template<class Type>
            void addSupType
            (
                const volScalarField& rho,
                const VolField<Type>& field,
                fvMatrix<Type>& eqn
            ) const;

            //- Add a source term to a phase equation
            template<class Type>
            void addSupType
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const VolField<Type>& field,
                fvMatrix<Type>& eqn
            ) const;


public:

    //- Runtime type information
    TypeName("fvTotalSource");


    // Constructors

        //- Construct from explicit source name and mesh
        fvTotalSource
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fvTotalSource();


    // Member Functions

        // Checks

            //- Return true if the fvModel adds a source term to the given
            //  field's transport equation
            virtual bool addsSupToField(const word& fieldName) const;


        // Access

            //- Return the phase name
            inline const word& phaseName() const;

            //- Return the volume of cells that the source applies to
            virtual scalar V() const = 0;


        // Sources

            //- Return the source value
            virtual dimensionedScalar S() const = 0;

            //- Return the source value
            virtual tmp<scalarField> source(const word& fieldName) const;


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvTotalSourceI.H"

#ifdef NoRepository
    #include "fvTotalSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
