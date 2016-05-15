/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::functionObjects::partialWrite

Group
    grpUtilitiesFunctionObjects

Description
    This function object allows user-selected fields/registered objects to be
    written at a custom write interval. The interval is given in terms of
    number of overall dumps

    Example of function object specification:
    \verbatim
    partialWrite1
    {
        type        partialWrite;
        functionObjectLibs ("libIOFunctionObjects.so");
        ...
        objectNames (p U T);
        writeInterval 100;
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: partialWrite | yes         |
        objectNames  | objects to write        | yes         |
        writeInterval | write interval         | yes         |
    \endtable

SeeAlso
    Foam::functionObject
    Foam::functionObjects::timeControl

SourceFiles
    partialWrite.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_partialWrite_H
#define functionObjects_partialWrite_H

#include "functionObject.H"
#include "HashSet.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class partialWrite Declaration
\*---------------------------------------------------------------------------*/

class partialWrite
:
    public functionObject
{
protected:

    // Protected data

        //- Reference to the objectRegistry
        const objectRegistry& obr_;

        //- Loaded fields
        UPtrList<volScalarField> vsf_;
        UPtrList<volVectorField> vvf_;
        UPtrList<volSphericalTensorField> vSpheretf_;
        UPtrList<volSymmTensorField> vSymmtf_;
        UPtrList<volTensorField> vtf_;

        UPtrList<surfaceScalarField> ssf_;
        UPtrList<surfaceVectorField> svf_;
        UPtrList<surfaceSphericalTensorField> sSpheretf_;
        UPtrList<surfaceSymmTensorField> sSymmtf_;
        UPtrList<surfaceTensorField> stf_;


        // Read from dictionary

            //- Names of objects to dump always
            HashSet<word> objectNames_;

            //- Write interval for restart dump
            label writeInterval_;


        //- Current dump instance. If reaches writeInterval do a full write.
        label writeInstance_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        partialWrite(const partialWrite&);

        //- Disallow default bitwise assignment
        void operator=(const partialWrite&);


        //- Load objects in the objectNames
        template<class Type>
        void loadField
        (
            const word&,
            UPtrList<GeometricField<Type, fvPatchField, volMesh>>&,
            UPtrList<GeometricField<Type, fvsPatchField, surfaceMesh>>&
        ) const;

        template<class Type>
        void changeWriteOptions
        (
            UPtrList<GeometricField<Type, fvPatchField, volMesh>>&,
            UPtrList<GeometricField<Type, fvsPatchField, surfaceMesh>>&,
            const IOobject::writeOption
        ) const;


public:

    //- Runtime type information
    TypeName("partialWrite");


    // Constructors

        //- Construct from Time and dictionary
        partialWrite
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~partialWrite();


    // Member Functions

        //- Read the partialWrite data
        virtual bool read(const dictionary&);

        //- Called when time was set at the end of the Time::operator++
        virtual bool timeSet();

        //- Execute
        virtual bool execute(const bool postProcess = false);

        //- Do nothing.
        //  The fields are registered and written automatically
        virtual bool write(const bool postProcess = false);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "partialWriteTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
