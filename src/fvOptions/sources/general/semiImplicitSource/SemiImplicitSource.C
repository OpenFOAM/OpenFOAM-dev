/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "SemiImplicitSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "DimensionedField.H"
#include "fvmSup.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template<class Type>
const Foam::wordList Foam::fv::SemiImplicitSource<Type>::volumeModeTypeNames_
(
    IStringStream("(absolute specific)")()
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
typename Foam::fv::SemiImplicitSource<Type>::volumeModeType
Foam::fv::SemiImplicitSource<Type>::wordToVolumeModeType
(
    const word& vmtName
) const
{
    forAll(volumeModeTypeNames_, i)
    {
        if (vmtName == volumeModeTypeNames_[i])
        {
            return volumeModeType(i);
        }
    }

    FatalErrorIn
    (
        "SemiImplicitSource<Type>::volumeModeType"
        "SemiImplicitSource<Type>::wordToVolumeModeType(const word&)"
    )   << "Unknown volumeMode type " << vmtName
        << ". Valid volumeMode types are:" << nl << volumeModeTypeNames_
        << exit(FatalError);

    return volumeModeType(0);
}


template<class Type>
Foam::word Foam::fv::SemiImplicitSource<Type>::volumeModeTypeToWord
(
    const volumeModeType& vmtType
) const
{
    if (vmtType > volumeModeTypeNames_.size())
    {
        return "UNKNOWN";
    }
    else
    {
        return volumeModeTypeNames_[vmtType];
    }
}


template<class Type>
void Foam::fv::SemiImplicitSource<Type>::setFieldData(const dictionary& dict)
{
    fieldNames_.setSize(dict.toc().size());
    injectionRate_.setSize(fieldNames_.size());

    applied_.setSize(fieldNames_.size(), false);

    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        fieldNames_[i] = iter().keyword();
        dict.lookup(iter().keyword()) >> injectionRate_[i];
        i++;
    }

    // Set volume normalisation
    if (volumeMode_ == vmAbsolute)
    {
        VDash_ = V_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::SemiImplicitSource<Type>::SemiImplicitSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    volumeMode_(vmAbsolute),
    VDash_(1.0),
    injectionRate_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fv::SemiImplicitSource<Type>::addSup
(
    fvMatrix<Type>& eqn,
    const label fieldI
)
{
    if (debug)
    {
        Info<< "SemiImplicitSource<" << pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    const GeometricField<Type, fvPatchField, volMesh>& psi = eqn.psi();

    DimensionedField<Type, volMesh> Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldI] + "Su",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<Type>
        (
            "zero",
            eqn.dimensions()/dimVolume,
            pTraits<Type>::zero
        ),
        false
    );

    UIndirectList<Type>(Su, cells_) = injectionRate_[fieldI].first()/VDash_;

    DimensionedField<scalar, volMesh> Sp
    (
        IOobject
        (
            name_ + fieldNames_[fieldI] + "Sp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<scalar>
        (
            "zero",
            Su.dimensions()/psi.dimensions(),
            0.0
        ),
        false
    );

    UIndirectList<scalar>(Sp, cells_) = injectionRate_[fieldI].second()/VDash_;

    eqn += Su + fvm::SuSp(Sp, psi);
}


template<class Type>
void Foam::fv::SemiImplicitSource<Type>::addSup
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const label fieldI
)
{
    if (debug)
    {
        Info<< "SemiImplicitSource<" << pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    return this->addSup(eqn, fieldI);
}


// ************************************************************************* //
