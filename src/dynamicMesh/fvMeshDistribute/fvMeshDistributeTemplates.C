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

        forAll(fld.boundaryField(), patchI)
        {
            Pout<< "    " << patchI
                << ' ' << fld.boundaryField()[patchI].patch().name()
                << ' ' << fld.boundaryField()[patchI].type()
                << ' ' << fld.boundaryField()[patchI].size()
                << endl;
        }
    }
}


// Save whole boundary field
template<class T, class Mesh>
void Foam::fvMeshDistribute::saveBoundaryFields
(
    PtrList<FieldField<fvsPatchField, T> >& bflds
) const
{
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


// Map boundary field
template<class T, class Mesh>
void Foam::fvMeshDistribute::mapBoundaryFields
(
    const mapPolyMesh& map,
    const PtrList<FieldField<fvsPatchField, T> >& oldBflds
)
{
    const labelList& oldPatchStarts = map.oldPatchStarts();
    const labelList& faceMap = map.faceMap();

    typedef GeometricField<T, fvsPatchField, Mesh> fldType;

    HashTable<fldType*> flds
    (
        mesh_.objectRegistry::lookupClass<fldType>()
    );

    if (flds.size() != oldBflds.size())
    {
        FatalErrorIn("fvMeshDistribute::mapBoundaryFields(..)") << "problem"
            << abort(FatalError);
    }

    label fieldI = 0;

    forAllIter(typename HashTable<fldType*>, flds, iter)
    {
        fldType& fld = *iter();
        typename fldType::GeometricBoundaryField& bfld =
            fld.boundaryField();

        const FieldField<fvsPatchField, T>& oldBfld = oldBflds[fieldI++];

        // Pull from old boundary field into bfld.

        forAll(bfld, patchI)
        {
            fvsPatchField<T>& patchFld = bfld[patchI];
            label faceI = patchFld.patch().start();

            forAll(patchFld, i)
            {
                label oldFaceI = faceMap[faceI++];

                // Find patch and local patch face oldFaceI was in.
                forAll(oldPatchStarts, oldPatchI)
                {
                    label oldLocalI = oldFaceI - oldPatchStarts[oldPatchI];

                    if (oldLocalI >= 0 && oldLocalI < oldBfld[oldPatchI].size())
                    {
                        patchFld[i] = oldBfld[oldPatchI][oldLocalI];
                    }
                }
            }
        }
    }
}


// Init patch fields of certain type
template<class GeoField, class PatchFieldType>
void Foam::fvMeshDistribute::initPatchFields
(
    const typename GeoField::value_type& initVal
)
{
    HashTable<GeoField*> flds
    (
        mesh_.objectRegistry::lookupClass<GeoField>()
    );

    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        GeoField& fld = *iter();

        typename GeoField::GeometricBoundaryField& bfld =
            fld.boundaryField();

        forAll(bfld, patchI)
        {
            if (isA<PatchFieldType>(bfld[patchI]))
            {
                bfld[patchI] == initVal;
            }
        }
    }
}


// correctBoundaryConditions patch fields of certain type
template<class GeoField>
void Foam::fvMeshDistribute::correctBoundaryConditions()
{
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


// Send fields. Note order supplied so we can receive in exactly the same order.
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
//
template<class GeoField>
void Foam::fvMeshDistribute::sendFields
(
    const label domain,
    const wordList& fieldNames,
    const fvMeshSubset& subsetter,
    Ostream& toNbr
)
{
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


// Opposite of sendFields
template<class GeoField>
void Foam::fvMeshDistribute::receiveFields
(
    const label domain,
    const wordList& fieldNames,
    fvMesh& mesh,
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
                    mesh.time().timeName(),
                    mesh,
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
