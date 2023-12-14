/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    const UPtrList<GeoField> fields(mesh.fields<GeoField>());

    forAll(fields, i)
    {
        const GeoField& field = fields[i];

        Pout<< "Field:" << field.name() << " internal size:" << field.size()
            << endl;

        forAll(field.boundaryField(), patchi)
        {
            Pout<< "    " << patchi
                << ' ' << field.boundaryField()[patchi].patch().name()
                << ' ' << field.boundaryField()[patchi].type()
                << ' ' << field.boundaryField()[patchi].size()
                << endl;
        }
    }
}


template<class Type, class Mesh>
void Foam::fvMeshDistribute::saveBoundaryFields
(
    PtrList<FieldField<fvsPatchField, Type>>& bfields
) const
{
    // Save whole boundary field

    const UPtrList<SurfaceField<Type>> fields
    (
        mesh_.fields<SurfaceField<Type>>()
    );

    bfields.setSize(fields.size());

    forAll(fields, i)
    {
        const SurfaceField<Type>& field = fields[i];

        bfields.set(i, field.boundaryField().clone().ptr());
    }
}


template<class Type, class Mesh>
void Foam::fvMeshDistribute::mapBoundaryFields
(
    const polyTopoChangeMap& map,
    const PtrList<FieldField<fvsPatchField, Type>>& oldBfields
)
{
    // Map boundary field

    const labelList& oldPatchStarts = map.oldPatchStarts();
    const labelList& faceMap = map.faceMap();

    UPtrList<SurfaceField<Type>> fields
    (
        mesh_.fields<SurfaceField<Type>>()
    );

    forAll(fields, i)
    {
        SurfaceField<Type>& field = fields[i];
        typename SurfaceField<Type>::Boundary& bfield =
            field.boundaryFieldRef();

        const FieldField<fvsPatchField, Type>& oldBfield = oldBfields[i];

        // Pull from old boundary field into bfield.

        forAll(bfield, patchi)
        {
            fvsPatchField<Type>& patchField = bfield[patchi];
            label facei = patchField.patch().start();

            forAll(patchField, i)
            {
                label oldFacei = faceMap[facei++];

                // Find patch and local patch face oldFacei was in.
                forAll(oldPatchStarts, oldPatchi)
                {
                    label oldLocalI = oldFacei - oldPatchStarts[oldPatchi];

                    if
                    (
                        oldLocalI >= 0
                     && oldLocalI < oldBfield[oldPatchi].size()
                    )
                    {
                        patchField[i] = oldBfield[oldPatchi][oldLocalI];
                    }
                }
            }
        }
    }
}


template<class Type>
void Foam::fvMeshDistribute::initMapExposedFaces
(
    PtrList<Field<Type>>& ifields
) const
{
    const UPtrList<SurfaceField<Type>> fields
    (
        mesh_.fields<SurfaceField<Type>>()
    );

    ifields.setSize(fields.size());

    forAll(fields, i)
    {
        ifields.set(i, Field<Type>(mesh_.nFaces()));

        const SurfaceField<Type>& field = fields[i];

        SubList<Type>(ifields[i], field.primitiveField().size()) =
            field.primitiveField();

        forAll(field.boundaryField(), patchi)
        {
            const fvsPatchField<Type>& pfield = field.boundaryField()[patchi];

            SubList<Type>(ifields[i], pfield.size(), pfield.patch().start()) =
                pfield;
        }
    }
}


template<class Type>
void Foam::fvMeshDistribute::mapExposedFaces
(
    const polyTopoChangeMap& map,
    const PtrList<Field<Type>>& oldFields
)
{
    UPtrList<SurfaceField<Type>> fields
    (
        mesh_.fields<SurfaceField<Type>>()
    );

    forAll(fields, i)
    {
        SurfaceField<Type>& field = fields[i];

        const Field<Type>& oldField = oldFields[i];

        const bool negateIfFlipped = isFlux(field);

        forAll(field.boundaryField(), patchi)
        {
            fvsPatchField<Type>& patchField = field.boundaryFieldRef()[patchi];

            forAll(patchField, i)
            {
                const label facei = patchField.patch().start()+i;
                const label oldFacei = map.faceMap()[facei];

                if (oldFacei < map.nOldInternalFaces())
                {
                    if (negateIfFlipped && map.flipFaceFlux().found(facei))
                    {
                        patchField[i] = flipOp()(oldField[oldFacei]);
                    }
                    else
                    {
                        patchField[i] = oldField[oldFacei];
                    }
                }
                else
                {
                    patchField[i] = oldField[oldFacei];
                }
            }
        }
    }
}


template<class GeoField>
void Foam::fvMeshDistribute::correctCoupledPatchFields()
{
    UPtrList<GeoField> fields(mesh_.fields<GeoField>());

    // Ensure the deltaCoeffs are available for constraint patch evaluation
    mesh_.deltaCoeffs();

    forAll(fields, i)
    {
        GeoField& field = fields[i];

        typename GeoField::Boundary& bfield = field.boundaryFieldRef();

        if
        (
            Pstream::defaultCommsType == Pstream::commsTypes::blocking
         || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            label nReq = Pstream::nRequests();

            forAll(bfield, patchi)
            {
                if (bfield[patchi].coupled())
                {
                    bfield[patchi].initEvaluate(Pstream::defaultCommsType);
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

            forAll(bfield, patchi)
            {
                if (bfield[patchi].coupled())
                {
                    bfield[patchi].evaluate(Pstream::defaultCommsType);
                }
            }
        }
        else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
        {
            const lduSchedule& patchSchedule =
                mesh_.globalData().patchSchedule();

            forAll(patchSchedule, patchEvali)
            {
                if (bfield[patchEvali].coupled())
                {
                    if (patchSchedule[patchEvali].init)
                    {
                        bfield[patchSchedule[patchEvali].patch]
                            .initEvaluate(Pstream::commsTypes::scheduled);
                    }
                    else
                    {
                        bfield[patchSchedule[patchEvali].patch]
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
        const GeoField& field =
            subsetter.baseMesh().lookupObject<GeoField>(fieldNames[i]);

        tmp<GeoField> tsubfield = subsetter.interpolate(field);

        toNbr
            << fieldNames[i] << token::NL << token::BEGIN_BLOCK
            << tsubfield
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
