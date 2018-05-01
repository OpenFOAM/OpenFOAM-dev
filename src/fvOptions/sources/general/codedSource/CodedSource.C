/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "CodedSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fv::CodedSource<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    word sourceType(pTraits<Type>::typeName);

    // Set additional rewrite rules
    dynCode.setFilterVariable("typeName", name_);
    dynCode.setFilterVariable("TemplateType", sourceType);
    dynCode.setFilterVariable("SourceType", sourceType + "Source");

    // dynCode.removeFilterVariable("code");
    dynCode.setFilterVariable("codeCorrect", codeCorrect_);
    dynCode.setFilterVariable("codeAddSup", codeAddSup_);
    dynCode.setFilterVariable("codeSetValue", codeSetValue_);

    // compile filtered C template
    dynCode.addCompileFile("codedFvOptionTemplate.C");

    // copy filtered H template
    dynCode.addCopyFile("codedFvOptionTemplate.H");

    // debugging: make BC verbose
    //         dynCode.setFilterVariable("verbose", "true");
    //         Info<<"compile " << name_ << " sha1: "
    //             << context.sha1() << endl;

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


template<class Type>
Foam::dlLibraryTable& Foam::fv::CodedSource<Type>::libs() const
{
    return const_cast<Time&>(mesh_.time()).libs();
}


template<class Type>
Foam::string Foam::fv::CodedSource<Type>::description() const
{
    return "fvOption:: " + name_;
}


template<class Type>
void Foam::fv::CodedSource<Type>::clearRedirect() const
{
    redirectFvOptionPtr_.clear();
}


template<class Type>
const Foam::dictionary& Foam::fv::CodedSource<Type>::codeDict() const
{
    return coeffs_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::CodedSource<Type>::CodedSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::fv::option& Foam::fv::CodedSource<Type>::redirectFvOption() const
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
void Foam::fv::CodedSource<Type>::correct
(
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    if (debug)
    {
        Info<< "CodedSource<"<< pTraits<Type>::typeName
            << ">::correct for source " << name_ << endl;
    }

    updateLibrary(name_);
    redirectFvOption().correct(field);
}


template<class Type>
void Foam::fv::CodedSource<Type>::addSup
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< "CodedSource<"<< pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    updateLibrary(name_);
    redirectFvOption().addSup(eqn, fieldi);
}


template<class Type>
void Foam::fv::CodedSource<Type>::addSup
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< "CodedSource<"<< pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    updateLibrary(name_);
    redirectFvOption().addSup(rho, eqn, fieldi);
}


template<class Type>
void Foam::fv::CodedSource<Type>::constrain
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< "CodedSource<"<< pTraits<Type>::typeName
            << ">::constrain for source " << name_ << endl;
    }

    updateLibrary(name_);
    redirectFvOption().constrain(eqn, fieldi);
}


// ************************************************************************* //
