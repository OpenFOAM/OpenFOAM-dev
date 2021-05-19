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

#include "codedFvModel.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(codedFvModel, 0);
    addToRunTimeSelectionTable(fvModel, codedFvModel, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::codedFvModel::readCoeffs()
{
    fieldName_ = coeffs().lookup<word>("field");

    if (fieldPrimitiveTypeName() != word::null)
    {
        updateLibrary();
    }
}


Foam::word Foam::fv::codedFvModel::fieldPrimitiveTypeName() const
{
    #define fieldPrimitiveTypeNameTernary(Type, nullArg)                       \
        mesh().foundObject<VolField<Type>>(fieldName_)                         \
      ? pTraits<Type>::typeName                                                \
      :

    return FOR_ALL_FIELD_TYPES(fieldPrimitiveTypeNameTernary) word::null;
}


void Foam::fv::codedFvModel::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    const word primitiveTypeName = fieldPrimitiveTypeName();

    // Set additional rewrite rules
    dynCode.setFilterVariable("typeName", name());
    dynCode.setFilterVariable("TemplateType", primitiveTypeName);
    dynCode.setFilterVariable("SourceType", primitiveTypeName + "Source");

    // compile filtered C template
    dynCode.addCompileFile("codedFvModelTemplate.C");

    // copy filtered H template
    dynCode.addCopyFile("codedFvModelTemplate.H");

    // define Make/options
    dynCode.setMakeOptions
    (
        "EXE_INC = -g \\\n"
        "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
        "-I$(LIB_SRC)/meshTools/lnInclude \\\n"
        "-I$(LIB_SRC)/sampling/lnInclude \\\n"
        "-I$(LIB_SRC)/fvModels/lnInclude \\\n"
        + context.options()
        + "\n\nLIB_LIBS = \\\n"
        + "    -lmeshTools \\\n"
        + "    -lfvModels \\\n"
        + "    -lsampling \\\n"
        + "    -lfiniteVolume \\\n"
        + context.libs()
    );
}


const Foam::word& Foam::fv::codedFvModel::codeName() const
{
    return name();
}


Foam::string Foam::fv::codedFvModel::description() const
{
    return "fvModel:: " + name();
}


void Foam::fv::codedFvModel::clearRedirect() const
{
    redirectFvModelPtr_.clear();
}


const Foam::dictionary& Foam::fv::codedFvModel::codeDict() const
{
    return coeffs();
}


Foam::wordList Foam::fv::codedFvModel::codeKeys() const
{

    return
    {
        "codeAddSup",
        "codeAddRhoSup",
        "codeAddAlphaRhoSup",
        "codeInclude",
        "localCode"
    };
}


Foam::fvModel& Foam::fv::codedFvModel::redirectFvModel() const
{
    if (!redirectFvModelPtr_.valid())
    {
        dictionary constructDict(coeffs());
        constructDict.set("type", name());
        redirectFvModelPtr_ = fvModel::New
        (
            name(),
            constructDict,
            mesh()
        );
    }
    return redirectFvModelPtr_();
}


template<class Type>
void Foam::fv::codedFvModel::addSupType
(
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    if (fieldPrimitiveTypeName() != word::null)
    {
        if (debug)
        {
            Info<< "codedFvModel::addSup for source " << name() << endl;
        }

        updateLibrary();
        redirectFvModel().addSup(eqn, fieldName);
    }
}


template<class Type>
void Foam::fv::codedFvModel::addSupType
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
            Info<< "codedFvModel::addSup for source " << name() << endl;
        }

        updateLibrary();
        redirectFvModel().addSup(rho, eqn, fieldName);
    }
}


template<class Type>
void Foam::fv::codedFvModel::addSupType
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
            Info<< "codedFvModel::addSup for source " << name() << endl;
        }

        updateLibrary();
        redirectFvModel().addSup(alpha, rho, eqn, fieldName);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::codedFvModel::codedFvModel
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    set_(coeffs(), mesh),
    fieldName_(word::null)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::codedFvModel::addSupFields() const
{
    return wordList(1, fieldName_);
}


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_SUP, fv::codedFvModel);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_RHO_SUP, fv::codedFvModel);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_SUP, fv::codedFvModel);


void Foam::fv::codedFvModel::updateMesh(const mapPolyMesh& mpm)
{
    set_.updateMesh(mpm);
}


bool Foam::fv::codedFvModel::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        set_.read(coeffs());
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
