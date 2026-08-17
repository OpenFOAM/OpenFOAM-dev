/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::LESfilter

Description
    Abstract class for LES filters

SourceFiles
    LESfilter.C
    newFilter.C

\*---------------------------------------------------------------------------*/

#ifndef LESfilter_H
#define LESfilter_H

#include "volFields.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "extrapolatedCalculatedFvPatchField.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                          Class LESfilter Declaration
\*---------------------------------------------------------------------------*/

class LESfilter
{
    // Private Data

        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("LESfilter");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            LESfilter,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& LESfilterDict
            ),
            (mesh, LESfilterDict)
        );


    // Constructors

        //- Construct from components
        LESfilter(const fvMesh& mesh)
        :
            mesh_(mesh)
        {}

        //- Disallow default bitwise copy construction
        LESfilter(const LESfilter&) = delete;


    // Selectors

        //- Return a reference to the selected LES filter
        static autoPtr<LESfilter> New
        (
            const fvMesh&,
            const dictionary&,
            const word& filterDictName="filter"
        );


    //- Destructor
    virtual ~LESfilter()
    {}


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Read the LESfilter dictionary
        virtual void read(const dictionary&) = 0;


    // Member Operators

        void operator=(const LESfilter&) = delete;

        virtual tmp<volInternalScalarField> operator[]
        (
            const volScalarField&
        ) const = 0;

        virtual tmp<volInternalVectorField> operator[]
        (
            const volVectorField&
        ) const = 0;

        virtual tmp<volInternalSymmTensorField> operator[]
        (
            const volSymmTensorField&
        ) const = 0;

        virtual tmp<volInternalTensorField> operator[]
        (
            const volTensorField&
        ) const = 0;

        template<class FieldType>
        tmp<typename FieldType::Internal> operator[]
        (
            const tmp<FieldType>& tf
        ) const
        {
            tmp<typename FieldType::Internal> result =  operator[](tf());
            tf.clear();
            return result;
        }

        template<class FieldType>
        tmp<FieldType> operator()(const FieldType& f) const
        {
            return FieldType::New
            (
                "filtered(" + f.name() + ')',
                operator[](f),
                extrapolatedCalculatedFvPatchField
                <
                    typename FieldType::Type_
                >::typeName
            );
        }

        template<class FieldType>
        tmp<FieldType> operator()(const tmp<FieldType>& tf) const
        {
            tmp<FieldType> result = FieldType::New
            (
                "filtered(" + tf().name() + ')',
                operator[](tf()),
                extrapolatedCalculatedFvPatchField
                <
                    typename FieldType::Type_
                >::typeName
            );

            tf.clear();

            return result;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
