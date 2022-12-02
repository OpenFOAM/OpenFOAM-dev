/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "polyTopoChangeMap.H"
#include "processorFvPatchField.H"

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

        Pout<< "Field:" << iter.key() << " internal size:" << fld.size()
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


template<class Type, class Mesh>
void Foam::fvMeshDistribute::saveBoundaryFields
(
    PtrList<FieldField<fvsPatchField, Type>>& bflds
) const
{
    // Save whole boundary field

    HashTable<const SurfaceField<Type>*> flds
    (
        static_cast<const fvMesh&>(mesh_)
       .objectRegistry::lookupClass<SurfaceField<Type>>()
    );

    bflds.setSize(flds.size());

    label i = 0;

    forAllConstIter(typename HashTable<const SurfaceField<Type>*>, flds, iter)
    {
        const SurfaceField<Type>& fld = *iter();

        bflds.set(i, fld.boundaryField().clone().ptr());

        i++;
    }
}


template<class Type, class Mesh>
void Foam::fvMeshDistribute::mapBoundaryFields
(
    const polyTopoChangeMap& map,
    const PtrList<FieldField<fvsPatchField, Type>>& oldBflds
)
{
    // Map boundary field

    const labelList& oldPatchStarts = map.oldPatchStarts();
    const labelList& faceMap = map.faceMap();

    HashTable<SurfaceField<Type>*> flds
    (
        mesh_.objectRegistry::lookupClass<SurfaceField<Type>>()
    );

    if (flds.size() != oldBflds.size())
    {
        FatalErrorInFunction
            << abort(FatalError);
    }

    label fieldi = 0;

    forAllIter(typename HashTable<SurfaceField<Type>*>, flds, iter)
    {
        SurfaceField<Type>& fld = *iter();
        typename SurfaceField<Type>::Boundary& bfld =
            fld.boundaryFieldRef();

        const FieldField<fvsPatchField, Type>& oldBfld = oldBflds[fieldi++];

        // Pull from old boundary field into bfld.

        forAll(bfld, patchi)
        {
            fvsPatchField<Type>& patchFld = bfld[patchi];
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


template<class Type>
void Foam::fvMeshDistribute::initMapExposedFaces
(
    PtrList<Field<Type>>& iflds
) const
{
    HashTable<const SurfaceField<Type>*> flds
    (
        static_cast<const fvMesh&>(mesh_).lookupClass<SurfaceField<Type>>()
    );

    iflds.setSize(flds.size());

    label fieldi = 0;

    forAllConstIter(typename HashTable<const SurfaceField<Type>*>, flds, iter)
    {
        iflds.set(fieldi, Field<Type>(mesh_.nFaces()));

        const SurfaceField<Type>& fld = *iter();

        SubList<Type>(iflds[fieldi], fld.primitiveField().size()) =
            fld.primitiveField();

        forAll(fld.boundaryField(), patchi)
        {
            const fvsPatchField<Type>& pfld = fld.boundaryField()[patchi];

            SubList<Type>(iflds[fieldi], pfld.size(), pfld.patch().start()) =
                pfld;
        }

        fieldi++;
    }
}


template<class Type>
void Foam::fvMeshDistribute::mapExposedFaces
(
    const polyTopoChangeMap& map,
    const PtrList<Field<Type>>& oldFlds
)
{
    HashTable<SurfaceField<Type>*> flds
    (
        mesh_.objectRegistry::lookupClass<SurfaceField<Type>>()
    );

    label fieldi = 0;

    forAllIter(typename HashTable<SurfaceField<Type>*>, flds, iter)
    {
        SurfaceField<Type>& fld = *iter();

        const Field<Type>& oldFld = oldFlds[fieldi];

        const bool negateIfFlipped = isFlux(fld);

        forAll(fld.boundaryField(), patchi)
        {
            fvsPatchField<Type>& patchFld = fld.boundaryFieldRef()[patchi];

            forAll(patchFld, i)
            {
                const label facei = patchFld.patch().start()+i;
                const label oldFacei = map.faceMap()[facei];

                if (oldFacei < map.nOldInternalFaces())
                {
                    if (negateIfFlipped && map.flipFaceFlux().found(facei))
                    {
                        patchFld[i] = flipOp()(oldFld[oldFacei]);
                    }
                    else
                    {
                        patchFld[i] = oldFld[oldFacei];
                    }
                }
                else
                {
                    patchFld[i] = oldFld[oldFacei];
                }
            }
        }

        fieldi++;
    }
}


template<class GeoField>
void Foam::fvMeshDistribute::correctCoupledPatchFields()
{
    HashTable<GeoField*> flds
    (
        mesh_.objectRegistry::lookupClass<GeoField>()
    );

    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        GeoField& fld = *iter();

        typename GeoField::Boundary& bfld = fld.boundaryFieldRef();

        if
        (
            Pstream::defaultCommsType == Pstream::commsTypes::blocking
         || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            label nReq = Pstream::nRequests();

            forAll(bfld, patchi)
            {
                if (bfld[patchi].coupled())
                {
                    bfld[patchi].initEvaluate(Pstream::defaultCommsType);
                }
            }

            // Block for any outstanding requests
            if
            (
                Pstream::parRun()
             && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
            )
            {
                Pstream::waitRequests(nReq);
            }

            forAll(bfld, patchi)
            {
                if (bfld[patchi].coupled())
                {
                    bfld[patchi].evaluate(Pstream::defaultCommsType);
                }
            }
        }
        else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
        {
            const lduSchedule& patchSchedule =
                mesh_.globalData().patchSchedule();

            forAll(patchSchedule, patchEvali)
            {
                if (bfld[patchEvali].coupled())
                {
                    if (patchSchedule[patchEvali].init)
                    {
                        bfld[patchSchedule[patchEvali].patch]
                            .initEvaluate(Pstream::commsTypes::scheduled);
                    }
                    else
                    {
                        bfld[patchSchedule[patchEvali].patch]
                            .evaluate(Pstream::commsTypes::scheduled);
                    }
                }
            }
        }
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
                    mesh.thisDb().time().name(),
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
