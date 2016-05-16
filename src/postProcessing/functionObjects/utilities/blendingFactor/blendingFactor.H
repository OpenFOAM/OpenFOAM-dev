/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::functionObjects::blendingFactor

Group
    grpUtilitiesFunctionObjects

Description
    This function object calculates and outputs the blendingFactor as used by
    the bended convection schemes.  The output is a volume field (cells) whose
    value is calculated via the maximum blending factor for any cell face.

SeeAlso
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    blendingFactor.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_blendingFactor_H
#define functionObjects_blendingFactor_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class blendingFactor Declaration
\*---------------------------------------------------------------------------*/

class blendingFactor
:
    public fvMeshFunctionObject
{
    // Private data

        //- Name of flux field, default is "phi"
        word phiName_;

        //- Field name
        word fieldName_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        blendingFactor(const blendingFactor&);

        //- Disallow default bitwise assignment
        void operator=(const blendingFactor&);

        //- Return the blending factor field from the database
        template<class Type>
        volScalarField& factor
        (
            const GeometricField<Type, fvPatchField, volMesh>& field
        );

        //- Calculate the blending factor
        template<class Type>
        void calc();


public:

    //- Runtime type information
    TypeName("blendingFactor");


    // Constructors

        //- Construct from Time and dictionary
        blendingFactor
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~blendingFactor();


    // Member Functions

        //- Read the blendingFactor data
        virtual bool read(const dictionary&);

        //- Calculate the blending-factor
        virtual bool execute(const bool postProcess = false);

        //- Write the blending-factor
        virtual bool write(const bool postProcess = false);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "blendingFactorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
