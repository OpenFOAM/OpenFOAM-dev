/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::readFields

Group
    grpFieldFunctionObjects

Description
    This function object reads fields from the time directories and adds them to
    the mesh database for further post-processing.

    Example of function object specification:
    \verbatim
    readFields1
    {
        type        readFields;
        functionObjectLibs ("libfieldFunctionObjects.so");
        ...
        fields
        (
            U
            p
        );
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: readFields   | yes         |
        fields       | list of fields to read  |  no         |
    \endtable

SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject

SourceFiles
    readFields.C
    IOreadFields.H

\*---------------------------------------------------------------------------*/

#ifndef readFields_H
#define readFields_H

#include "OFstream.H"
#include "pointFieldFwd.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                         Class readFields Declaration
\*---------------------------------------------------------------------------*/

class readFields
{
protected:

    // Protected data

        //- Name of this set of readFields object
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Fields to load
        wordList fieldSet_;

        //- Loaded fields
        PtrList<volScalarField> vsf_;
        PtrList<volVectorField> vvf_;
        PtrList<volSphericalTensorField> vSpheretf_;
        PtrList<volSymmTensorField> vSymmtf_;
        PtrList<volTensorField> vtf_;

        PtrList<surfaceScalarField> ssf_;
        PtrList<surfaceVectorField> svf_;
        PtrList<surfaceSphericalTensorField> sSpheretf_;
        PtrList<surfaceSymmTensorField> sSymmtf_;
        PtrList<surfaceTensorField> stf_;


    // Protected Member Functions

        //- Disallow default bitwise copy construct
        readFields(const readFields&);

        //- Disallow default bitwise assignment
        void operator=(const readFields&);

        template<class Type>
        void loadField
        (
            const word&,
            PtrList<GeometricField<Type, fvPatchField, volMesh> >&,
            PtrList<GeometricField<Type, fvsPatchField, surfaceMesh> >&
        ) const;


public:

    //- Runtime type information
    TypeName("readFields");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        readFields
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~readFields();


    // Member Functions

        //- Return name of the readFields object
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the field min/max data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "readFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
