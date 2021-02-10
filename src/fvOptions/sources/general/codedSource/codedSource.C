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

\*---------------------------------------------------------------------------*/

#include "codedSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Private Static Data Members * * * * * * * * * * * //

const Foam::wordList Foam::fv::codedSource::codeKeys_ =
{
    "codeAddSup",
    "codeAddRhoSup",
    "codeAddAlphaRhoSup",
    "codeInclude",
    "localCode"
};


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(codedSource, 0);
    addToRunTimeSelectionTable(option, codedSource, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::codedSource::readCoeffs()
{
    fieldName_ = coeffs_.lookup<word>("field");

    name_ = coeffs_.lookup<word>("name");

    if (fieldPrimitiveTypeName() != word::null)
    {
        updateLibrary();
    }
}


Foam::word Foam::fv::codedSource::fieldPrimitiveTypeName() const
{
    #define fieldPrimitiveTypeNameTernary(Type, nullArg)                       \
        mesh_.foundObject<VolField<Type>>(fieldName_)                          \
      ? pTraits<Type>::typeName                                                \
      :

    return FOR_ALL_FIELD_TYPES(fieldPrimitiveTypeNameTernary) word::null;
}


void Foam::fv::codedSource::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    const word primitiveTypeName = fieldPrimitiveTypeName();

    // Set additional rewrite rules
    dynCode.setFilterVariable("typeName", name_);
    dynCode.setFilterVariable("TemplateType", primitiveTypeName);
    dynCode.setFilterVariable("SourceType", primitiveTypeName + "Source");

    // compile filtered C template
    dynCode.addCompileFile("codedFvOptionTemplate.C");

    // copy filtered H template
    dynCode.addCopyFile("codedFvOptionTemplate.H");

    // define Make/options
    dynCode.setMakeOptions
    (
        "EXE_INC = -g \\\n"
        "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
        "-I$(LIB_SRC)/meshTools/lnInclude \\\n"
        "-I$(LIB_SRC)/sampling/lnInclude \\\n"
        "-I$(LIB_SRC)/fvOptions/lnInclude \\\n"
        + context.options()
        + "\n\nLIB_LIBS = \\\n"
        + "    -lmeshTools \\\n"
        + "    -lfvOptions \\\n"
        + "    -lsampling \\\n"
        + "    -lfiniteVolume \\\n"
        + context.libs()
    );
}


const Foam::word& Foam::fv::codedSource::codeName() const
{
    return name_;
}


Foam::string Foam::fv::codedSource::description() const
{
    return "fvOption:: " + name_;
}


void Foam::fv::codedSource::clearRedirect() const
{
    redirectFvOptionPtr_.clear();
}


const Foam::dictionary& Foam::fv::codedSource::codeDict() const
{
    return coeffs_;
}


const Foam::wordList& Foam::fv::codedSource::codeKeys() const
{
    return codeKeys_;
}


Foam::fv::option& Foam::fv::codedSource::redirectFvOption() const
{
    if (!redirectFvOptionPtr_.valid())
    {
        dictionary constructDict(dict_);
        constructDict.set("type", name_);
        redirectFvOptionPtr_ = option::New
        (
            name_,
            constructDict,
            mesh_
        );
    }
    return redirectFvOptionPtr_();
}


template<class Type>
void Foam::fv::codedSource::addSupType
(
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    if (fieldPrimitiveTypeName() != word::null)
    {
        if (debug)
        {
            Info<< "codedSource::addSup for source " << name_ << endl;
        }

        updateLibrary();
        redirectFvOption().addSup(eqn, fieldName);
    }
}


template<class Type>
void Foam::fv::codedSource::addSupType
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    if (fieldPrimitiveTypeName() != word::null)
    {
        if (debug)
        {
            Info<< "codedSource::addSup for source " << name_ << endl;
        }

        updateLibrary();
        redirectFvOption().addSup(rho, eqn, fieldName);
    }
}


template<class Type>
void Foam::fv::codedSource::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    if (fieldPrimitiveTypeName() != word::null)
    {
        if (debug)
        {
            Info<< "codedSource::addSup for source " << name_ << endl;
        }

        updateLibrary();
        redirectFvOption().addSup(alpha, rho, eqn, fieldName);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::codedSource::codedSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    fieldName_(word::null)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::codedSource::addSupFields() const
{
    return wordList(1, fieldName_);
}


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_OPTION_ADD_SUP, codedSource);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_OPTION_ADD_RHO_SUP, codedSource);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_OPTION_ADD_ALPHA_RHO_SUP, codedSource);


bool Foam::fv::codedSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
