/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::functionObjects::fieldMinMax

Description
    Calculates the value and location of scalar minimum and maximum for a list
    of user-specified fields.

    For variables with a rank greater than zero, either the min/max of a
    component value or the magnitude is reported.  When operating in parallel,
    the processor owning the value is also given.

    Example of function object specification:
    \verbatim
    fieldMinMax1
    {
        type        fieldMinMax;
        libs        ("libfieldFunctionObjects.so");
        ...
        log         yes;
        location    yes;
        mode        magnitude;
        fields
        (
            U
            p
        );
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: fieldMinMax  | yes         |
        log          | write min/max data to standard output | no | no
        location     | write location of the min/max value | no | yes
        mode         | calculation mode: magnitude or component | no | magnitude
    \endtable

    Output data is written to the file \<timeDir\>/fieldMinMax.dat

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles

SourceFiles
    fieldMinMax.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldMinMax_H
#define functionObjects_fieldMinMax_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class fieldMinMax Declaration
\*---------------------------------------------------------------------------*/

class fieldMinMax
:
    public fvMeshFunctionObject,
    public logFiles
{
public:

    enum class modeType
    {
        mag,
        cmpt
    };

protected:

    // Protected data

        //- Mode type names
        static const NamedEnum<modeType, 2> modeTypeNames_;

        //- Switch to write location of min/max values
        Switch location_;

        //- Mode for min/max - only applicable for ranks > 0
        modeType mode_;

        //- Fields to assess min/max
        wordList fieldSet_;


    // Protected Member Functions

        //- Helper function to write the output
        template<class Type>
        void output
        (
            const word& fieldName,
            const word& outputName,
            const label minCell,
            const label maxCell,
            const vector& minC,
            const vector& maxC,
            const label minProci,
            const label maxProci,
            const Type& minValue,
            const Type& maxValue
        );

        //- Calculate the field min/max
        template<class Type>
        void calcMinMaxFields
        (
            const word& fieldName,
            const modeType& mode
        );

        //- Output file header information
        virtual void writeFileHeader(const label i);


public:

    //- Runtime type information
    TypeName("fieldMinMax");


    // Constructors

        //- Construct from Time and dictionary
        fieldMinMax
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        fieldMinMax(const fieldMinMax&) = delete;


    //- Destructor
    virtual ~fieldMinMax();


    // Member Functions

        //- Read the field min/max data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the fieldMinMax
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fieldMinMax&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldMinMaxTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
