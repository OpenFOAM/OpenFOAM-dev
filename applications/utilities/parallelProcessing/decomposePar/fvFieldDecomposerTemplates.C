/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "fvFieldDecomposer.H"
#include "processorFvPatchField.H"
#include "processorFvsPatchField.H"
#include "processorCyclicFvPatchField.H"
#include "processorCyclicFvsPatchField.H"
#include "emptyFvPatchFields.H"
#include "stringOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvFieldDecomposer::mapCellToFace
(
    const labelUList& owner,
    const labelUList& neighbour,
    const Field<Type>& field,
    const labelUList& addressing
)
{
    tmp<Field<Type>> tfld(new Field<Type>(addressing.size()));
    Field<Type>& fld = tfld.ref();

    forAll(addressing, i)
    {
        fld[i] =
            field
            [
                addressing[i] > 0
              ? neighbour[addressing[i] - 1]
              : owner[- addressing[i] - 1]
            ];
    }

    return tfld;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvFieldDecomposer::mapFaceToFace
(
    const Field<Type>& field,
    const labelUList& addressing,
    const bool isFlux
)
{
    tmp<Field<Type>> tfld(new Field<Type>(addressing.size()));
    Field<Type>& fld = tfld.ref();

    forAll(addressing, i)
    {
        fld[i] =
            (isFlux && addressing[i] < 0 ? -1 : +1)
           *field[mag(addressing[i]) - 1];
    }

    return tfld;
}


template<class Type>
Foam::PtrList<typename Foam::VolField<Type>::Internal>
Foam::fvFieldDecomposer::decomposeVolInternalField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field
    const typename VolField<Type>::Internal field
    (
        IOobject
        (
            fieldIoObject.name(),
            completeMesh_.time().name(),
            completeMesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        completeMesh_
    );

    // Construct the processor fields
    PtrList<typename VolField<Type>::Internal> procFields(procMeshes_.size());
    forAll(procMeshes_, proci)
    {
        // Create the processor field with the dummy patch fields
        procFields.set
        (
            proci,
            new typename VolField<Type>::Internal
            (
                IOobject
                (
                    fieldIoObject.name(),
                    procMeshes_[proci].time().name(),
                    procMeshes_[proci],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                procMeshes_[proci],
                field.dimensions(),
                Field<Type>(field.primitiveField(), cellProcAddressing_[proci])
            )
        );
    }

    return procFields;
}


template<class Type>
Foam::PtrList<Foam::VolField<Type>>
Foam::fvFieldDecomposer::decomposeVolField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field
    const VolField<Type> field
    (
        IOobject
        (
            fieldIoObject.name(),
            completeMesh_.time().name(),
            completeMesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        completeMesh_
    );

    // Construct the processor fields
    PtrList<VolField<Type>> procFields(procMeshes_.size());
    forAll(procMeshes_, proci)
    {
        // Create dummy patch fields
        PtrList<fvPatchField<Type>> patchFields
        (
            procMeshes_[proci].boundary().size()
        );
        forAll(procMeshes_[proci].boundary(), procPatchi)
        {
            patchFields.set
            (
                procPatchi,
                fvPatchField<Type>::New
                (
                    calculatedFvPatchField<Type>::typeName,
                    procMeshes_[proci].boundary()[procPatchi],
                    DimensionedField<Type, volMesh>::null()
                )
            );
        }

        // Create the processor field with the dummy patch fields
        procFields.set
        (
            proci,
            new VolField<Type>
            (
                IOobject
                (
                    fieldIoObject.name(),
                    procMeshes_[proci].time().name(),
                    procMeshes_[proci],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                procMeshes_[proci],
                field.dimensions(),
                Field<Type>(field.primitiveField(), cellProcAddressing_[proci]),
                patchFields,
                field.sources().table()
            )
        );

        // Alias the created proc field
        VolField<Type>& vf = procFields[proci];

        // Change the patch fields to the correct type using a mapper
        // constructor (with reference to the now correct internal field)
        typename VolField<Type>::Boundary& bf = vf.boundaryFieldRef();
        forAll(bf, procPatchi)
        {
            const fvPatch& procPatch =
                procMeshes_[proci].boundary()[procPatchi];

            const label completePatchi = completePatchID(proci, procPatchi);

            if (completePatchi == procPatchi)
            {
                bf.set
                (
                    procPatchi,
                    fvPatchField<Type>::New
                    (
                        field.boundaryField()[completePatchi],
                        procPatch,
                        vf(),
                        patchFieldDecomposers_[proci][procPatchi]
                    )
                );
            }
            else if (isA<processorCyclicFvPatch>(procPatch))
            {
                if (field.boundaryField()[completePatchi].overridesConstraint())
                {
                    OStringStream str;
                    str << "\nThe field \"" << field.name()
                        << "\" on cyclic patch \""
                        << field.boundaryField()[completePatchi].patch().name()
                        << "\" cannot be decomposed as it is not a cyclic "
                        << "patch field. A \"patchType cyclic;\" setting has "
                        << "been used to override the cyclic patch type.\n\n"
                        << "Cyclic patches like this with non-cyclic boundary "
                        << "conditions should be confined to a single "
                        << "processor using decomposition constraints.";
                    FatalErrorInFunction
                        << stringOps::breakIntoIndentedLines(str.str()).c_str()
                        << exit(FatalError);
                }

                const label nbrCompletePatchi =
                    refCast<const processorCyclicFvPatch>(procPatch)
                   .referPatch().nbrPatchIndex();

                // Use `fvPatchField<Type>::New` rather than
                // `new processorCyclicFvPatchField<Type>` so that derivations
                // (such as non-conformal processor cyclics) are constructed
                bf.set
                (
                    procPatchi,
                    fvPatchField<Type>::New
                    (
                        procPatch.type(),
                        procPatch,
                        vf()
                    )
                );

                bf[procPatchi] =
                    mapCellToFace
                    (
                        labelUList(),
                        completeMesh_.lduAddr().patchAddr(nbrCompletePatchi),
                        field.primitiveField(),
                        faceProcAddressingBf_[proci][procPatchi]
                    );
            }
            else if (isA<processorFvPatch>(procPatch))
            {
                bf.set
                (
                    procPatchi,
                    fvPatchField<Type>::New
                    (
                        procPatch.type(),
                        procPatch,
                        vf()
                    )
                );

                bf[procPatchi] =
                    mapCellToFace
                    (
                        completeMesh_.owner(),
                        completeMesh_.neighbour(),
                        field.primitiveField(),
                        faceProcAddressingBf_[proci][procPatchi]
                    );
            }
            else
            {
                FatalErrorInFunction
                    << "Unknown type." << abort(FatalError);
            }
        }
    }

    return procFields;
}


template<class Type>
Foam::PtrList<Foam::SurfaceField<Type>>
Foam::fvFieldDecomposer::decomposeFvSurfaceField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field
    const SurfaceField<Type> field
    (
        IOobject
        (
            fieldIoObject.name(),
            completeMesh_.time().name(),
            completeMesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        completeMesh_
    );

    // Construct the processor fields
    PtrList<SurfaceField<Type>> procFields(procMeshes_.size());
    forAll(procMeshes_, proci)
    {
        const SubList<label> faceAddressingIf
        (
            faceProcAddressing_[proci],
            procMeshes_[proci].nInternalFaces()
        );

        // Create dummy patch fields
        PtrList<fvsPatchField<Type>> patchFields
        (
            procMeshes_[proci].boundary().size()
        );
        forAll(procMeshes_[proci].boundary(), procPatchi)
        {
            patchFields.set
            (
                procPatchi,
                fvsPatchField<Type>::New
                (
                    calculatedFvsPatchField<Type>::typeName,
                    procMeshes_[proci].boundary()[procPatchi],
                    DimensionedField<Type, surfaceMesh>::null()
                )
            );
        }

        // Create the processor field with the dummy patch fields
        procFields.set
        (
            proci,
            new SurfaceField<Type>
            (
                IOobject
                (
                    field.name(),
                    procMeshes_[proci].time().name(),
                    procMeshes_[proci],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                procMeshes_[proci],
                field.dimensions(),
                mapFaceToFace
                (
                    field,
                    faceAddressingIf,
                    isFlux(field)
                ),
                patchFields
            )
        );

        // Alias the created proc field
        SurfaceField<Type>& sf = procFields[proci];

        // Change the patch fields to the correct type using a mapper
        // constructor (with reference to the now correct internal field)
        typename SurfaceField<Type>::Boundary& bf = sf.boundaryFieldRef();
        forAll(procMeshes_[proci].boundary(), procPatchi)
        {
            const fvPatch& procPatch =
                procMeshes_[proci].boundary()[procPatchi];

            const label completePatchi = completePatchID(proci, procPatchi);

            if (completePatchi == procPatchi)
            {
                bf.set
                (
                    procPatchi,
                    fvsPatchField<Type>::New
                    (
                        field.boundaryField()[procPatchi],
                        procPatch,
                        sf(),
                        patchFieldDecomposers_[proci][procPatchi]
                    )
                );
            }
            else if (isA<processorCyclicFvPatch>(procPatch))
            {
                bf.set
                (
                    procPatchi,
                    new processorCyclicFvsPatchField<Type>
                    (
                        procPatch,
                        sf(),
                        mapFaceToFace
                        (
                            field.boundaryField()[completePatchi],
                            faceProcAddressingBf_[proci][procPatchi],
                            isFlux(field)
                        )
                    )
                );
            }
            else if (isA<processorFvPatch>(procPatch))
            {
                bf.set
                (
                    procPatchi,
                    new processorFvsPatchField<Type>
                    (
                        procPatch,
                        sf(),
                        mapFaceToFace
                        (
                            field.primitiveField(),
                            faceProcAddressingBf_[proci][procPatchi],
                            isFlux(field)
                        )
                    )
                );
            }
            else
            {
                FatalErrorInFunction
                    << "Unknown type." << abort(FatalError);
            }
        }
    }

    return procFields;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvFieldDecomposer::decomposeVolInternalFields
(
    const IOobjectList& objects
)
{
    const word& fieldClassName = VolField<Type>::Internal::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< nl << "    Decomposing " << fieldClassName << "s" << nl << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            Info<< "        " << fieldIter()->name() << endl;

            PtrList<typename VolField<Type>::Internal> procFields =
                decomposeVolInternalField<Type>(*fieldIter());

            forAll(procFields, proci)
            {
                procFields[proci].write();
            }
        }
    }
}


template<class Type>
void Foam::fvFieldDecomposer::decomposeVolFields
(
    const IOobjectList& objects
)
{
    const word& fieldClassName = VolField<Type>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< nl << "    Decomposing " << fieldClassName << "s" << nl << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            Info<< "        " << fieldIter()->name() << endl;

            PtrList<VolField<Type>> procFields =
                decomposeVolField<Type>(*fieldIter());

            forAll(procFields, proci)
            {
                procFields[proci].write();
            }
        }
    }
}


template<class Type>
void Foam::fvFieldDecomposer::decomposeFvSurfaceFields
(
    const IOobjectList& objects
)
{
    const word& fieldClassName = SurfaceField<Type>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< nl << "    Decomposing " << fieldClassName << "s" << nl << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            Info<< "        " << fieldIter()->name() << endl;

            PtrList<SurfaceField<Type>> procFields =
                decomposeFvSurfaceField<Type>(*fieldIter());

            forAll(procFields, proci)
            {
                procFields[proci].write();
            }
        }
    }
}


// ************************************************************************* //
