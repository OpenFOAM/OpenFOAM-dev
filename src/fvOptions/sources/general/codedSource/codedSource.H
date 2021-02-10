/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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
    Foam::fv::codedSource

Description
   Constructs on-the-fly fvOption source from user-supplied code

Usage
    Example usage in constant/fvOptions:
    \verbatim
    energySource
    {
        type            codedSource;
        name            sourceTime;

        selectionMode   all;

        field           h;

        codeInclude
        #{
        #};

        codeAddSup
        #{
            Pout<< "**codeAddSup**" << endl;
            const Time& time = mesh().time();
            const scalarField& V = mesh_.V();
            scalarField& heSource = eqn.source();
            heSource -= 0.1*sqr(time.value())*V;
        #};

        codeAddRhoSup
        #{
            Pout<< "**codeAddRhoSup**" << endl;
        #};

        codeAddAlphaRhoSup
        #{
            Pout<< "**codeAddAlphaRhoSup**" << endl;
        #};
    }
    \endverbatim


SourceFiles
    codedSource.C

\*---------------------------------------------------------------------------*/

#ifndef codedSource_H
#define codedSource_H

#include "cellSetOption.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class codedSource Declaration
\*---------------------------------------------------------------------------*/

class codedSource
:
    public cellSetOption,
    public codedBase
{
    // Private static data

        //- The keywords associated with source code
        static const wordList codeKeys_;


    // Private data

        //- The name
        word name_;

        //- The name of the field that this option applies to
        word fieldName_;

        //- Underlying functionObject
        mutable autoPtr<option> redirectFvOptionPtr_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- Return the name of the field's primitive type
        word fieldPrimitiveTypeName() const;

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;

        //- Name of the dynamically generated CodedType
        virtual const word& codeName() const;

        //- Return a description (type + name) for the output
        virtual string description() const;

        //- Clear any redirected objects
        virtual void clearRedirect() const;

        //- Get the dictionary to initialize the codeContext
        virtual const dictionary& codeDict() const;

        //- Get the keywords associated with source code
        virtual const wordList& codeKeys() const;

        //- Dynamically compiled fvOption
        option& redirectFvOption() const;


        // Sources

            //- Add a source term to an equation
            template<class Type>
            void addSupType
            (
                fvMatrix<Type>& eqn,
                const word& fieldName
            ) const;

            //- Add a source term to a compressible equation
            template<class Type>
            void addSupType
            (
                const volScalarField& rho,
                fvMatrix<Type>& eqn,
                const word& fieldName
            ) const;

            //- Add a source term to a phase equation
            template<class Type>
            void addSupType
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<Type>& eqn,
                const word& fieldName
            ) const;


public:

    //- Runtime type information
    TypeName("codedSource");


    // Constructors

        //- Construct from components
        codedSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Sources

            //- Add a source term to an equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_OPTION_ADD_SUP);

            //- Add a source term to a compressible equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_OPTION_ADD_RHO_SUP);

            //- Add a source term to a phase equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_OPTION_ADD_ALPHA_RHO_SUP);


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
