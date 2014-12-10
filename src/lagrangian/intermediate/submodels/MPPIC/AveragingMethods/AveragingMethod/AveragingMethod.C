/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "AveragingMethod.H"
#include "runTimeSelectionTables.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::AveragingMethod<Type>::AveragingMethod
(
    const IOobject& io,
    const dictionary& dict,
    const fvMesh& mesh,
    const labelList& size
)
:
    regIOobject(io),
    FieldField<Field, Type>(),
    dict_(dict),
    mesh_(mesh)
{
    forAll(size, i)
    {
        FieldField<Field, Type>::append
        (
            new Field<Type>(size[i], pTraits<Type>::zero)
        );
    }
}


template<class Type>
Foam::AveragingMethod<Type>::AveragingMethod
(
    const AveragingMethod<Type>& am
)
:
    regIOobject(am),
    FieldField<Field, Type>(am),
    dict_(am.dict_),
    mesh_(am.mesh_)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::AveragingMethod<Type> >
Foam::AveragingMethod<Type>::New
(
    const IOobject& io,
    const dictionary& dict,
    const fvMesh& mesh
)
{
    word averageType(dict.lookup(typeName));

    //Info<< "Selecting averaging method "
    //    << averageType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(averageType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "Foam::AveragingMethod<Type>::New"
            "("
                "const IOobject&, "
                "const dictionary&, "
                "const fvMesh&"
            ")"
        )   << "Unknown averaging method " << averageType
            << ", constructor not in hash table" << nl << nl
            << "    Valid averaging methods are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<AveragingMethod<Type> >(cstrIter()(io, dict, mesh));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::AveragingMethod<Type>::~AveragingMethod()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::AveragingMethod<Type>::updateGrad()
{
    // do nothing
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::AveragingMethod<Type>::average()
{
    updateGrad();
}


template<class Type>
void Foam::AveragingMethod<Type>::average
(
    const AveragingMethod<scalar>& weight
)
{
    updateGrad();

    *this /= max(weight, SMALL);
}


template<class Type>
bool Foam::AveragingMethod<Type>::writeData(Ostream& os) const
{
    return os.good();
}


template<class Type>
bool Foam::AveragingMethod<Type>::write() const
{
    const pointMesh pointMesh_(mesh_);

    // point volumes
    Field<scalar> pointVolume(mesh_.nPoints(), 0);

    // output fields
    GeometricField<Type, fvPatchField, volMesh> cellValue
    (
        IOobject
        (
            this->name() + ":cellValue",
            this->time().timeName(),
            mesh_
        ),
        mesh_,
        dimensioned<Type>("zero", dimless, pTraits<Type>::zero)
    );
    GeometricField<TypeGrad, fvPatchField, volMesh> cellGrad
    (
        IOobject
        (
            this->name() + ":cellGrad",
            this->time().timeName(),
            mesh_
        ),
        mesh_,
        dimensioned<TypeGrad>("zero", dimless, pTraits<TypeGrad>::zero)
    );
    GeometricField<Type, pointPatchField, pointMesh> pointValue
    (
        IOobject
        (
            this->name() + ":pointValue",
            this->time().timeName(),
            mesh_
        ),
        pointMesh_,
        dimensioned<Type>("zero", dimless, pTraits<Type>::zero)
    );
    GeometricField<TypeGrad, pointPatchField, pointMesh> pointGrad
    (
        IOobject
        (
            this->name() + ":pointGrad",
            this->time().timeName(),
            mesh_
        ),
        pointMesh_,
        dimensioned<TypeGrad>("zero", dimless, pTraits<TypeGrad>::zero)
    );

    // tet-volume weighted sums
    forAll(mesh_.C(), cellI)
    {
        const List<tetIndices> cellTets =
            polyMeshTetDecomposition::cellTetIndices(mesh_, cellI);

        forAll(cellTets, tetI)
        {
            const tetIndices& tetIs = cellTets[tetI];
            const scalar v = tetIs.tet(mesh_).mag();

            cellValue[cellI] += v*interpolate(mesh_.C()[cellI], tetIs);
            cellGrad[cellI] += v*interpolateGrad(mesh_.C()[cellI], tetIs);

            const face& f = mesh_.faces()[tetIs.face()];
            labelList vertices(3);
            vertices[0] = f[tetIs.faceBasePt()];
            vertices[1] = f[tetIs.facePtA()];
            vertices[2] = f[tetIs.facePtB()];

            forAll(vertices, vertexI)
            {
                const label pointI = vertices[vertexI];

                pointVolume[pointI] += v;
                pointValue[pointI] +=
                    v*interpolate(mesh_.points()[pointI], tetIs);
                pointGrad[pointI] +=
                    v*interpolateGrad(mesh_.points()[pointI], tetIs);
            }
        }
    }

    // average
    cellValue.internalField() /= mesh_.V();
    cellGrad.internalField() /= mesh_.V();
    pointValue.internalField() /= pointVolume;
    pointGrad.internalField() /= pointVolume;

    // write
    if (!cellValue.write()) return false;
    if (!cellGrad.write()) return false;
    if (!pointValue.write()) return false;
    if (!pointGrad.write()) return false;

    return true;
}


// ************************************************************************* //
