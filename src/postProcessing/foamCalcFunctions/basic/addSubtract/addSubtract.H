/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::calcTypes::addSubtract

Description
    adds/subtracts a field or value to/from a base field.

    New field name specified by -resultName option, or automatically as:
        \<baseFieldName\>_add_<addSubtractFieldName>
        \<baseFieldName\>_add_value
        \<baseFieldName\>_subtract_<addSubtractFieldName>
        \<baseFieldName\>_subtract_value

    Example usage:
        addSubtract p add -value 100000 -resultName pAbs
        addSubtract U subtract -field U0

SourceFiles
    addSubtract.C
    writeaddSubtractField.C
    writeaddSubtractValue.C

\*---------------------------------------------------------------------------*/

#ifndef addSubtract_H
#define addSubtract_H

#include "calcType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace calcTypes
{

/*---------------------------------------------------------------------------*\
                          Class addSubtract Declaration
\*---------------------------------------------------------------------------*/

class addSubtract
:
    public calcType
{
public:

    enum calcTypes
    {
        FIELD,
        VALUE
    };

    enum calcModes
    {
        ADD,
        SUBTRACT
    };


private:

    // Private data

        //- Name of base field (to addSubtract to)
        word baseFieldName_;

        //- Calc type as given by enumerations above
        calcTypes calcType_;

        //- Name of field to add/subtract
        word addSubtractFieldName_;

        //- String representation of value to add/subtract
        string addSubtractValueStr_;

        //- Name of result field
        word resultName_;

        //- Mode - addSubtract/subtract
        calcModes calcMode_;


    // Private Member Functions

        // Output

            //- Calc and output field addSubtractitions
            void writeAddSubtractFields
            (
                const Time& runTime,
                const fvMesh& mesh,
                const IOobject& baseFieldHeader
            );

            //- Calc and output field and value addSubtractitions
            void writeAddSubtractValues
            (
                const Time& runTime,
                const fvMesh& mesh,
                const IOobject& baseFieldHeader
            );


        //- Disallow default bitwise copy construct
        addSubtract(const addSubtract&);

        //- Disallow default bitwise assignment
        void operator=(const addSubtract&);


protected:

    // Member Functions

        // Calculation routines

            //- Initialise - typically setting static variables,
            //  e.g. command line arguments
            virtual void init();

            //- Pre-time loop calculations
            virtual void preCalc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );

            //- Time loop calculations
            virtual void calc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );


        // I-O

            //- Write addSubtract field
            template<class Type>
            void writeAddSubtractField
            (
                const IOobject& baseHeader,
                const IOobject& addSubtractHeader,
                const fvMesh& mesh,
                bool& processed
            );

            //- Write addSubtract value
            template<class Type>
            void writeAddSubtractValue
            (
                const IOobject& baseHeader,
                const string& valueStr,
                const fvMesh& mesh,
                bool& processed
            );


public:

    //- Runtime type information
    TypeName("addSubtract");


    // Constructors

        //- Construct null
        addSubtract();


    //- Destructor
    virtual ~addSubtract();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace calcTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "writeAddSubtractField.C"
#   include "writeAddSubtractValue.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
