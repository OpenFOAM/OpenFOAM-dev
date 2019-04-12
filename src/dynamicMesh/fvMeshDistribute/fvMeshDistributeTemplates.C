/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "mapPolyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoField>
void Foam::fvMeshDistribute::printFieldInfo(const fvMesh& mesh)
{
    HashTable<const GeoField*> flds
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    forAllConstIter(typename HashTable<const GeoField*>, flds, iter)
    {
        const GeoField& fld = *iter();

        Pout<< "Field:" << iter.key() << " internalsize:" << fld.size()
            //<< " value:" << fld
            << endl;

        forAll(fld.boundaryField(), patchi)
        {
            Pout<< "    " << patchi
                << ' ' << fld.boundaryField()[patchi].patch().name()
                << ' ' << fld.boundaryField()[patchi].type()
                << ' ' << fld.boundaryField()[patchi].size()
                << endl;
        }
    }
}


template<class T, class Mesh>
void Foam::fvMeshDistribute::saveBoundaryFields
(
    PtrList<FieldField<fvsPatchField, T>>& bflds
) const
{
    // Save whole boundary field

    typedef GeometricField<T, fvsPatchField, Mesh> fldType;

    HashTable<const fldType*> flds
    (
        static_cast<const fvMesh&>(mesh_).objectRegistry::lookupClass<fldType>()
    );

    bflds.setSize(flds.size());

    label i = 0;

    forAllConstIter(typename HashTable<const fldType*>, flds, iter)
    {
        const fldType& fld = *iter();

        bflds.set(i, fld.boundaryField().clone().ptr());

        i++;
    }
}


template<class T, class Mesh>
void Foam::fvMeshDistribute::mapBoundaryFields
(
    const mapPolyMesh& map,
    const PtrList<FieldField<fvsPatchField, T>>& oldBflds
)
{
    // Map boundary field

    const labelList& oldPatchStarts = map.oldPatchStarts();
    const labelList& faceMap = map.faceMap();

    typedef GeometricField<T, fvsPatchField, Mesh> fldType;

    HashTable<fldType*> flds
    (
        mesh_.objectRegistry::lookupClass<fldType>()
    );

    if (flds.size() != oldBflds.size())
    {
        FatalErrorInFunction
            << abort(FatalError);
    }

    label fieldi = 0;

    forAllIter(typename HashTable<fldType*>, flds, iter)
    {
        fldType& fld = *iter();
        typename fldType::Boundary& bfld =
            fld.boundaryFieldRef();

        const FieldField<fvsPatchField, T>& oldBfld = oldBflds[fieldi++];

        // Pull from old boundary field into bfld.

        forAll(bfld, patchi)
        {
            fvsPatchField<T>& patchFld = bfld[patchi];
            label facei = patchFld.patch().start();

            forAll(patchFld, i)
            {
                label oldFacei = faceMap[facei++];

                // Find patch and local patch face oldFacei was in.
                forAll(oldPatchStarts, oldPatchi)
                {
                    label oldLocalI = oldFacei - oldPatchStarts[oldPatchi];

                    if (oldLocalI >= 0 && oldLocalI < oldBfld[oldPatchi].size())
                    {
                        patchFld[i] = oldBfld[oldPatchi][oldLocalI];
                    }
                }
            }
        }
    }
}


template<class T>
void Foam::fvMeshDistribute::saveInternalFields
(
    PtrList<Field<T>>& iflds
) const
{
    typedef GeometricField<T, fvsPatchField, surfaceMesh> fldType;

    HashTable<const fldType*> flds
    (
        static_cast<const fvMesh&>(mesh_).objectRegistry::lookupClass<fldType>()
    );

    iflds.setSize(flds.size());

    label i = 0;

    forAllConstIter(typename HashTable<const fldType*>, flds, iter)
    {
        const fldType& fld = *iter();

        iflds.set(i, fld.primitiveField().clone());

        i++;
    }
}


template<class T>
void Foam::fvMeshDistribute::mapExposedFaces
(
    const mapPolyMesh& map,
    const PtrList<Field<T>>& oldFlds
)
{
    // Set boundary values of exposed internal faces

    const labelList& faceMap = map.faceMap();

    typedef GeometricField<T, fvsPatchField, surfaceMesh> fldType;

    HashTable<fldType*> flds
    (
        mesh_.objectRegistry::lookupClass<fldType>()
    );

    if (flds.size() != oldFlds.size())
    {
        FatalErrorIn("fvMeshDistribute::mapExposedFaces(..)") << "problem"
            << abort(FatalError);
    }


    label fieldI = 0;

    forAllIter(typename HashTable<fldType*>, flds, iter)
    {
        fldType& fld = *iter();
        typename fldType::Boundary& bfld = fld.boundaryFieldRef();

        const Field<T>& oldInternal = oldFlds[fieldI++];

        // Pull from old internal field into bfld.

        forAll(bfld, patchi)
        {
            fvsPatchField<T>& patchFld = bfld[patchi];

            forAll(patchFld, i)
            {
                const label faceI = patchFld.patch().start()+i;

                label oldFaceI = faceMap[faceI];

                if (oldFaceI < oldInternal.size())
                {
                    patchFld[i] = oldInternal[oldFaceI];

                    if (map.flipFaceFlux().found(faceI))
                    {
                        patchFld[i] = flipOp()(patchFld[i]);
                    }
                }
            }
        }
    }
}


template<class GeoField, class PatchFieldType>
void Foam::fvMeshDistribute::initPatchFields
(
    const typename GeoField::value_type& initVal
)
{
    // Init patch fields of certain type

    HashTable<GeoField*> flds
    (
        mesh_.objectRegistry::lookupClass<GeoField>()
    );

    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        GeoField& fld = *iter();

        typename GeoField::Boundary& bfld = fld.boundaryFieldRef();

        forAll(bfld, patchi)
        {
            if (isA<PatchFieldType>(bfld[patchi]))
            {
                bfld[patchi] == initVal;
            }
        }
    }
}


template<class GeoField>
void Foam::fvMeshDistribute::correctBoundaryConditions()
{
    // correctBoundaryConditions patch fields of certain type

    HashTable<GeoField*> flds
    (
        mesh_.objectRegistry::lookupClass<GeoField>()
    );

    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        const GeoField& fld = *iter();
        fld.correctBoundaryConditions();
    }
}


template<class GeoField>
void Foam::fvMeshDistribute::sendFields
(
    const label domain,
    const wordList& fieldNames,
    const fvMeshSubset& subsetter,
    Ostream& toNbr
)
{
    // Send fields. Note order supplied so we can receive in exactly the same
    // order.
    // Note that field gets written as entry in dictionary so we
    // can construct from subdictionary.
    // (since otherwise the reading as-a-dictionary mixes up entries from
    // consecutive fields)
    // The dictionary constructed is:
    //  volScalarField
    //  {
    //      p {internalField ..; boundaryField ..;}
    //      k {internalField ..; boundaryField ..;}
    //  }
    //  volVectorField
    //  {
    //      U {internalField ...  }
    //  }

    // volVectorField {U {internalField ..; boundaryField ..;}}

    toNbr << GeoField::typeName << token::NL << token::BEGIN_BLOCK << token::NL;
    forAll(fieldNames, i)
    {
        if (debug)
        {
            Pout<< "Subsetting field " << fieldNames[i]
                << " for domain:" << domain << endl;
        }

        // Send all fieldNames. This has to be exactly the same set as is
        // being received!
        const GeoField& fld =
            subsetter.baseMesh().lookupObject<GeoField>(fieldNames[i]);

        tmp<GeoField> tsubfld = subsetter.interpolate(fld);

        toNbr
            << fieldNames[i] << token::NL << token::BEGIN_BLOCK
            << tsubfld
            << token::NL << token::END_BLOCK << token::NL;
    }
    toNbr << token::END_BLOCK << token::NL;
}


template<class GeoField>
void Foam::fvMeshDistribute::receiveFields
(
    const label domain,
    const wordList& fieldNames,
    typename GeoField::Mesh& mesh,
    PtrList<GeoField>& fields,
    const dictionary& fieldDicts
)
{
    if (debug)
    {
        Pout<< "Receiving fields " << fieldNames
            << " from domain:" << domain << endl;
    }

    fields.setSize(fieldNames.size());

    forAll(fieldNames, i)
    {
        if (debug)
        {
            Pout<< "Constructing field " << fieldNames[i]
                << " from domain:" << domain << endl;
        }

        fields.set
        (
            i,
            new GeoField
            (
                IOobject
                (
                    fieldNames[i],
                    mesh.thisDb().time().timeName(),
                    mesh.thisDb(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                fieldDicts.subDict(fieldNames[i])
            )
        );
    }
}


// ************************************************************************* //
