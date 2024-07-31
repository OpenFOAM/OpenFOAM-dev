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
    Foam::fvConstraint

Description
    Finite volume options abstract base class.

SourceFiles
    fvConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef fvConstraint_H
#define fvConstraint_H

#include "fvMatricesFwd.H"
#include "volFieldsFwd.H"
#include "dictionary.H"
#include "fvConstraintM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class polyTopoChangeMap;
class polyMeshMap;
class polyDistributionMap;

/*---------------------------------------------------------------------------*\
                           Class option Declaration
\*---------------------------------------------------------------------------*/

class fvConstraint
{
    // Private Member Data

        //- Constraint name
        const word name_;

        //- Constraint type
        const word constraintType_;

        //- Reference to the mesh database
        const fvMesh& mesh_;


    // Private Member Functions

        //- Apply a constraint to an equation
        template<class Type>
        bool constrainType
        (
            fvMatrix<Type>& eqn,
            const word& fieldName
        ) const;

        //- Apply constraint to a field
        template<class Type>
        bool constrainType(VolField<Type>& field) const;


public:

    //- Runtime type information
    TypeName("fvConstraint");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fvConstraint,
            dictionary,
            (
                const word& name,
                const word& constraintType,
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (name, constraintType, mesh, dict)
        );


    // Constructors

        //- Construct from components
        fvConstraint
        (
            const word& name,
            const word& constraintType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Return clone
        autoPtr<fvConstraint> clone() const
        {
            NotImplemented;
            return autoPtr<fvConstraint>(nullptr);
        }

        //- Return pointer to new fvConstraint object created
        //  on the freestore from an Istream
        class iNew
        {
            //- Reference to the mesh
            const fvMesh& mesh_;

            const word& name_;

        public:

            iNew
            (
                const fvMesh& mesh,
                const word& name
            )
            :
                mesh_(mesh),
                name_(name)
            {}

            autoPtr<fvConstraint> operator()(Istream& is) const
            {
                // const word name(is);
                const dictionary dict(is);

                return autoPtr<fvConstraint>
                (
                    fvConstraint::New(name_, mesh_, dict)
                );
            }
        };


    // Selectors

        //- Return a reference to the selected fvConstraint
        static autoPtr<fvConstraint> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fvConstraint();


    // Member Functions

        // Access

            //- Return const access to the source name
            inline const word& name() const;

            //- Return const access to the mesh database
            inline const fvMesh& mesh() const;

            //- Return the coefficients sub-dictionary
            inline const dictionary& coeffs(const dictionary&) const;


        // Checks

            //- Return the list of constrained fields
            virtual wordList constrainedFields() const;

            //- Return true if the given field is constrained
            virtual bool constrainsField(const word& fieldName) const;


        // Constraints

            //- Apply a constraint to an equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_CONSTRAINT_CONSTRAIN);

            //- Apply constraint to a field
            FOR_ALL_FIELD_TYPES(DEFINE_FV_CONSTRAINT_CONSTRAIN_FIELD);


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints() = 0;

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&) = 0;

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&) = 0;

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&) = 0;


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvConstraintI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
