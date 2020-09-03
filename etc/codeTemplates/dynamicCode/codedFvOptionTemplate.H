/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

Description
    Template for use with dynamic code generation of a source.

SourceFiles
    codedFvOptionTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef codedFvOptionTemplate_H
#define codedFvOptionTemplate_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                         A templated CodedFvOption
\*---------------------------------------------------------------------------*/

class ${typeName}FvOption${SourceType}
:
    public cellSetOption
{
public:

    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;

    //- Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from patch and internal field
        ${typeName}FvOption${SourceType}
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~${typeName}FvOption${SourceType}();


    // Member Functions

            //- Correct field
            virtual void correct
            (
                GeometricField<${TemplateType}, fvPatchField, volMesh>&
            ) const;

            //- Explicit and implicit matrix contributions
            virtual void addSup
            (
                fvMatrix<${TemplateType}>& eqn,
                const label fieldi
            ) const;

            //- Explicit and implicit matrix contributions for compressible
            //  equations
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<${TemplateType}>& eqn,
                const label fieldi
            ) const;

            //- Explicit and implicit matrix contributions for phase equations
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<${TemplateType}>& eqn,
                const label fieldi
            ) const;

            //- Set value
            virtual void constrain
            (
                fvMatrix<${TemplateType}>& eqn,
                const label fieldi
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
