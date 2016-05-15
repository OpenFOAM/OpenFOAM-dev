/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::functionObjects::grad

Group
    grpFieldFunctionObjects

Description
    This function object calculates the gradient of a field.  The operation is
    limited to scalar and vector volume or surface fields, and the output is a
    volume vector or tensor field.

SourceFiles
    grad.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_grad_H
#define functionObjects_grad_H

#include "functionObject.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dimensionSet;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class grad Declaration
\*---------------------------------------------------------------------------*/

class grad
:
    public functionObject
{
    // Private data

        //- Reference to the objectRegistry
        const objectRegistry& obr_;

        //- Name of field to process
        word fieldName_;

        //- Name of result field
        word resultName_;


    // Private Member Functions

        //- Helper function to create/store/return the gradient field
        template<class Type>
        GeometricField
        <
            typename outerProduct<vector, Type>::type,
            fvPatchField,
            volMesh
        >&
        gradField(const word& gradName, const dimensionSet& dims);

        //- Helper function to calculate the gradient of different field types
        template<class Type>
        void calcGrad
        (
            const word& fieldName,
            const word& resultName,
            bool& processed
        );

        //- Disallow default bitwise copy construct
        grad(const grad&);

        //- Disallow default bitwise assignment
        void operator=(const grad&);


public:

    //- Runtime type information
    TypeName("grad");


    // Constructors

        //- Construct from Time and dictionary
        grad
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~grad();


    // Member Functions

        //- Read the grad data
        virtual bool read(const dictionary&);

        //- Calculate the gradient field
        virtual bool execute(const bool postProcess = false);

        //- Write the gradient field
        virtual bool write(const bool postProcess = false);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "gradTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
