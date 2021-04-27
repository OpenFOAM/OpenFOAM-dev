/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::fvConstraints

Description
    Finite volume constraints

SourceFiles
    fvConstraints.C

\*---------------------------------------------------------------------------*/

#ifndef fvConstraints_H
#define fvConstraints_H

#include "fvConstraint.H"
#include "PtrListDictionary.H"
#include "MeshObject.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class fvConstraints;

Ostream& operator<<(Ostream& os, const fvConstraints& constraints);


/*---------------------------------------------------------------------------*\
                           Class fvConstraints Declaration
\*---------------------------------------------------------------------------*/

class fvConstraints
:
    public MeshObject<fvMesh, UpdateableMeshObject, fvConstraints>,
    public dictionary,
    public PtrListDictionary<fvConstraint>
{
    // Private Member Data

        //- Time index to check that all defined sources have been applied
        mutable label checkTimeIndex_;

        //- Sets of the fields that have been constrained by the fvConstraints
        mutable PtrList<wordHashSet> constrainedFields_;


    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- Check that all fvConstraints have been applied
        void checkApplied() const;


public:

    //- Runtime type information
    TypeName("fvConstraints");


    // Constructors

        //- Construct from components with list of field names
        fvConstraints(const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        fvConstraints(const fvConstraints&) = delete;

        //- Inherit the base New method
        using MeshObject<fvMesh, UpdateableMeshObject, fvConstraints>::New;


    //- Destructor
    virtual ~fvConstraints()
    {}


    // Member Functions

        //- Declare fvConstraints to be a global dictionary
        virtual bool global() const
        {
            return true;
        }

        //- Return complete path + object name if the file exists
        //  either in the case/processor or case otherwise null
        virtual fileName filePath() const
        {
            return globalFilePath(type());
        }


        // Checks

            //- Return true if an fvConstraint constrains the given field
            virtual bool constrainsField(const word& fieldName) const;


        // Constraints

            //- Apply constraints to an equation
            template<class Type>
            bool constrain(fvMatrix<Type>& eqn) const;

            //- Apply constraints to a field
            template<class Type>
            bool constrain
            (
                GeometricField<Type, fvPatchField, volMesh>& field
            ) const;


        // Mesh changes

            //- Update for mesh changes
            virtual void updateMesh(const mapPolyMesh&);

            //- Update for mesh motion
            virtual bool movePoints();


        // IO

            //- ReadData function which reads the fvConstraints dictionary
            //  required for regIOobject read operation
            virtual bool readData(Istream&);

            //- Write data to Ostream
            virtual bool writeData(Ostream& os) const;

            //- Read the fvConstraints dictionary if it has changed
            //  and update the constraints
            virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvConstraints&) = delete;


    // IOstream Operators

        //- Ostream operator
        friend Ostream& operator<<
        (
            Ostream& os,
            const fvConstraints& constraints
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Template function for obtaining global status
template<>
inline bool typeGlobal<fvConstraints>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvConstraintsTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
