/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "GeometricField.H"
#include "Time.H"
#include "demandDrivenData.H"
#include "dictionary.H"
#include "data.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// check mesh for two fields

#define checkField(gf1, gf2, op)                                    \
if ((gf1).mesh() != (gf2).mesh())                                   \
{                                                                   \
    FatalErrorIn("checkField(gf1, gf2, op)")                        \
        << "different mesh for fields "                             \
        << (gf1).name() << " and " << (gf2).name()                  \
        << " during operatrion " <<  op                             \
        << abort(FatalError);                                       \
}


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::readFields
(
    const dictionary& dict
)
{
    DimensionedField<Type, GeoMesh>::readField(dict, "internalField");

    boundaryField_.readField(*this, dict.subDict("boundaryField"));

    if (dict.found("referenceLevel"))
    {
        Type fieldAverage(pTraits<Type>(dict.lookup("referenceLevel")));

        Field<Type>::operator+=(fieldAverage);

        forAll(boundaryField_, patchi)
        {
            boundaryField_[patchi] == boundaryField_[patchi] + fieldAverage;
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::readFields()
{
    const IOdictionary dict
    (
        IOobject
        (
            this->name(),
            this->time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->readStream(typeName)
    );

    this->close();

    readFields(dict);
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::readIfPresent()
{
    if
    (
        this->readOpt() == IOobject::MUST_READ
     || this->readOpt() == IOobject::MUST_READ_IF_MODIFIED
    )
    {
        WarningIn
        (
            "GeometricField<Type, PatchField, GeoMesh>::readIfPresent()"
        )   << "read option IOobject::MUST_READ or MUST_READ_IF_MODIFIED"
            << " suggests that a read constructor for field " << this->name()
            << " would be more appropriate." << endl;
    }
    else if (this->readOpt() == IOobject::READ_IF_PRESENT && this->headerOk())
    {
        readFields();

        // Check compatibility between field and mesh
        if (this->size() != GeoMesh::size(this->mesh()))
        {
            FatalIOErrorIn
            (
                "GeometricField<Type, PatchField, GeoMesh>::"
                "readIfPresent()",
                this->readStream(typeName)
            )   << "   number of field elements = " << this->size()
                << " number of mesh elements = "
                << GeoMesh::size(this->mesh())
                << exit(FatalIOError);
        }

        readOldTimeIfPresent();

        return true;
    }

    return false;
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::readOldTimeIfPresent()
{
    // Read the old time field if present
    IOobject field0
    (
        this->name()  + "_0",
        this->time().timeName(),
        this->db(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE,
        this->registerObject()
    );

    if (field0.headerOk())
    {
        if (debug)
        {
            Info<< "Reading old time level for field"
                << endl << this->info() << endl;
        }

        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            field0,
            this->mesh()
        );

        field0Ptr_->timeIndex_ = timeIndex_ - 1;

        if (!field0Ptr_->readOldTimeIfPresent())
        {
            field0Ptr_->oldTime();
        }

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

// Constructor given a GeometricField and dimensionSet
// This allocates storage for the field but not values.
// Note : This constructor should only be used to
//       construct TEMPORARY variables
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const word& patchFieldType
)
:
    DimensionedField<Type, GeoMesh>(io, mesh, ds, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(mesh.boundary(), *this, patchFieldType)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "creating temporary"
            << endl << this->info() << endl;
    }

    readIfPresent();
}


// Constructor given a GeometricField and dimensionSet
// This allocates storage for the field but not values.
// Note : This constructor should only be used to
//       construct TEMPORARY variables
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
:
    DimensionedField<Type, GeoMesh>(io, mesh, ds, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(mesh.boundary(), *this, patchFieldTypes, actualPatchTypes)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "creating temporary"
            << endl << this->info() << endl;
    }

    readIfPresent();
}


// Constructor given a GeometricField and dimensioned<Type>
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const word& patchFieldType
)
:
    DimensionedField<Type, GeoMesh>(io, mesh, dt, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(mesh.boundary(), *this, patchFieldType)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "creating temporary"
            << endl << this->info() << endl;
    }

    boundaryField_ == dt.value();

    readIfPresent();
}


// Constructor given a GeometricField and dimensioned<Type>
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
:
    DimensionedField<Type, GeoMesh>(io, mesh, dt, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(mesh.boundary(), *this, patchFieldTypes, actualPatchTypes)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "creating temporary"
            << endl << this->info() << endl;
    }

    boundaryField_ == dt.value();

    readIfPresent();
}


//  construct from components
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const Field<Type>& iField,
    const PtrList<PatchField<Type> >& ptfl
)
:
    DimensionedField<Type, GeoMesh>(io, mesh, ds, iField),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(mesh.boundary(), *this, ptfl)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "constructing from components"
            << endl << this->info() << endl;
    }

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh
)
:
    DimensionedField<Type, GeoMesh>(io, mesh, dimless, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(mesh.boundary())
{
    readFields();

    // Check compatibility between field and mesh

    if (this->size() != GeoMesh::size(this->mesh()))
    {
        FatalIOErrorIn
        (
            "GeometricField<Type, PatchField, GeoMesh>::GeometricField"
            "(const IOobject&, const Mesh&)",
            this->readStream(typeName)
        )   << "   number of field elements = " << this->size()
            << " number of mesh elements = " << GeoMesh::size(this->mesh())
            << exit(FatalIOError);
    }

    readOldTimeIfPresent();

    if (debug)
    {
        Info<< "Finishing read-construct of "
               "GeometricField<Type, PatchField, GeoMesh>"
            << endl << this->info() << endl;
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dictionary& dict
)
:
    DimensionedField<Type, GeoMesh>(io, mesh, dimless, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(mesh.boundary())
{
    readFields(dict);

    // Check compatibility between field and mesh

    if (this->size() != GeoMesh::size(this->mesh()))
    {
        FatalErrorIn
        (
            "GeometricField<Type, PatchField, GeoMesh>::GeometricField"
            "(const IOobject&, const Mesh&, const dictionary&)"
        )   << "   number of field elements = " << this->size()
            << " number of mesh elements = " << GeoMesh::size(this->mesh())
            << exit(FatalIOError);
    }

    readOldTimeIfPresent();

    if (debug)
    {
        Info<< "Finishing dictionary-construct of "
               "GeometricField<Type, PatchField, GeoMesh>"
            << endl << this->info() << endl;
    }
}


// construct as copy
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
:
    DimensionedField<Type, GeoMesh>(gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(*this, gf.boundaryField_)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "constructing as copy"
            << endl << this->info() << endl;
    }

    if (gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            *gf.field0Ptr_
        );
    }

    this->writeOpt() = IOobject::NO_WRITE;
}

// construct as copy of tmp<GeometricField> deleting argument
#ifdef ConstructFromTmp
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf
)
:
    DimensionedField<Type, GeoMesh>
    (
        const_cast<GeometricField<Type, PatchField, GeoMesh>&>(tgf()),
        tgf.isTmp()
    ),
    timeIndex_(tgf().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(*this, tgf().boundaryField_)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "constructing as copy"
            << endl << this->info() << endl;
    }

    this->writeOpt() = IOobject::NO_WRITE;

    tgf.clear();
}
#endif


// construct as copy resetting IO parameters
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
:
    DimensionedField<Type, GeoMesh>(io, gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(*this, gf.boundaryField_)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "constructing as copy resetting IO params"
            << endl << this->info() << endl;
    }

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            io.name() + "_0",
            *gf.field0Ptr_
        );
    }
}


// construct as copy resetting name
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const word& newName,
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
:
    DimensionedField<Type, GeoMesh>(newName, gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(*this, gf.boundaryField_)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "constructing as copy resetting name"
            << endl << this->info() << endl;
    }

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            newName + "_0",
            *gf.field0Ptr_
        );
    }
}


// construct as copy resetting name
#ifdef ConstructFromTmp
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const word& newName,
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf
)
:
    DimensionedField<Type, GeoMesh>
    (
        newName,
        const_cast<GeometricField<Type, PatchField, GeoMesh>&>(tgf()),
        tgf.isTmp()
    ),
    timeIndex_(tgf().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(*this, tgf().boundaryField_)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "constructing from tmp resetting name"
            << endl << this->info() << endl;
    }

    tgf.clear();
}
#endif

// construct as copy resetting IO parameters and patch type
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, PatchField, GeoMesh>& gf,
    const word& patchFieldType
)
:
    DimensionedField<Type, GeoMesh>(io, gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_(this->mesh().boundary(), *this, patchFieldType)
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "constructing as copy resetting IO params"
            << endl << this->info() << endl;
    }

    boundaryField_ == gf.boundaryField_;

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            io.name() + "_0",
            *gf.field0Ptr_
        );
    }
}


// construct as copy resetting IO parameters and boundary types
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, PatchField, GeoMesh>& gf,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes

)
:
    DimensionedField<Type, GeoMesh>(io, gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    boundaryField_
    (
        this->mesh().boundary(),
        *this,
        patchFieldTypes,
        actualPatchTypes
    )
{
    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "constructing as copy resetting IO params"
            << endl << this->info() << endl;
    }

    boundaryField_ == gf.boundaryField_;

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            io.name() + "_0",
            *gf.field0Ptr_
        );
    }
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::~GeometricField()
{
    deleteDemandDrivenData(field0Ptr_);
    deleteDemandDrivenData(fieldPrevIterPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
typename
Foam::GeometricField<Type, PatchField, GeoMesh>::DimensionedInternalField&
Foam::GeometricField<Type, PatchField, GeoMesh>::dimensionedInternalField()
{
    this->setUpToDate();
    storeOldTimes();
    return *this;
}


template<class Type, template<class> class PatchField, class GeoMesh>
typename
Foam::GeometricField<Type, PatchField, GeoMesh>::InternalField&
Foam::GeometricField<Type, PatchField, GeoMesh>::internalField()
{
    this->setUpToDate();
    storeOldTimes();
    return *this;
}


// Return reference to GeometricBoundaryField
template<class Type, template<class> class PatchField, class GeoMesh>
typename
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField&
Foam::GeometricField<Type, PatchField, GeoMesh>::boundaryField()
{
    this->setUpToDate();
    storeOldTimes();
    return boundaryField_;
}


// Store old-time field
template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::storeOldTimes() const
{
    if
    (
        field0Ptr_
     && timeIndex_ != this->time().timeIndex()
     && !(
            this->name().size() > 2
         && this->name()(this->name().size()-2, 2) == "_0"
         )
    )
    {
        storeOldTime();
    }

    // Correct time index
    timeIndex_ = this->time().timeIndex();
}

// Store old-time field
template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::storeOldTime() const
{
    if (field0Ptr_)
    {
        field0Ptr_->storeOldTime();

        if (debug)
        {
            Info<< "Storing old time field for field" << endl
                << this->info() << endl;
        }

        *field0Ptr_ == *this;
        field0Ptr_->timeIndex_ = timeIndex_;

        if (field0Ptr_->field0Ptr_)
        {
            field0Ptr_->writeOpt() = this->writeOpt();
        }
    }
}

// Return the number of old time fields stored
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::label Foam::GeometricField<Type, PatchField, GeoMesh>::nOldTimes() const
{
    if (field0Ptr_)
    {
        return field0Ptr_->nOldTimes() + 1;
    }
    else
    {
        return 0;
    }
}

// Return old time internal field
template<class Type, template<class> class PatchField, class GeoMesh>
const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::oldTime() const
{
    if (!field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                this->name() + "_0",
                this->time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                this->registerObject()
            ),
            *this
        );
    }
    else
    {
        storeOldTimes();
    }

    return *field0Ptr_;
}

// Return old time internal field
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::oldTime()
{
    static_cast<const GeometricField<Type, PatchField, GeoMesh>&>(*this)
        .oldTime();

    return *field0Ptr_;
}


// Store previous iteration field
template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::storePrevIter() const
{
    if (!fieldPrevIterPtr_)
    {
        if (debug)
        {
            Info<< "Allocating previous iteration field" << endl
                << this->info() << endl;
        }

        fieldPrevIterPtr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            this->name() + "PrevIter",
            *this
        );
    }
    else
    {
        *fieldPrevIterPtr_ == *this;
    }
}


// Return previous iteration field
template<class Type, template<class> class PatchField, class GeoMesh>
const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::prevIter() const
{
    if (!fieldPrevIterPtr_)
    {
        FatalErrorIn
        (
            "GeometricField<Type, PatchField, GeoMesh>::prevIter() const"
        )   << "previous iteration field" << endl << this->info() << endl
            << "  not stored."
            << "  Use field.storePrevIter() at start of iteration."
            << abort(FatalError);
    }

    return *fieldPrevIterPtr_;
}


// Correct the boundary conditions
template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::
correctBoundaryConditions()
{
    this->setUpToDate();
    storeOldTimes();
    boundaryField_.evaluate();
}


// Does the field need a reference level for solution
template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::needReference() const
{
    // Search all boundary conditions, if any are
    // fixed-value or mixed (Robin) do not set reference level for solution.

    bool needRef = true;

    forAll(boundaryField_, patchi)
    {
        if (boundaryField_[patchi].fixesValue())
        {
            needRef = false;
            break;
        }
    }

    reduce(needRef, andOp<bool>());

    return needRef;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::relax(const scalar alpha)
{
    if (debug)
    {
        InfoIn
        (
            "GeometricField<Type, PatchField, GeoMesh>::relax"
            "(const scalar alpha)"
        )  << "Relaxing" << endl << this->info() << " by " << alpha << endl;
    }

    operator==(prevIter() + alpha*(*this - prevIter()));
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::relax()
{
    word name = this->name();

    if
    (
        this->mesh().data::template lookupOrDefault<bool>
        (
            "finalIteration",
            false
        )
    )
    {
        name += "Final";
    }

    if (this->mesh().relaxField(name))
    {
        relax(this->mesh().fieldRelaxationFactor(name));
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::word Foam::GeometricField<Type, PatchField, GeoMesh>::select
(
    bool final
) const
{
    if (final)
    {
        return this->name() + "Final";
    }
    else
    {
        return this->name();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::writeMinMax
(
    Ostream& os
) const
{
    os  << "min/max(" << this->name() << ") = "
        << Foam::min(*this).value() << ", "
        << Foam::max(*this).value()
        << endl;
}


// writeData member function required by regIOobject
template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::
writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh> >
Foam::GeometricField<Type, PatchField, GeoMesh>::T() const
{
    tmp<GeometricField<Type, PatchField, GeoMesh> > result
    (
        new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                this->name() + ".T()",
                this->instance(),
                this->db()
            ),
            this->mesh(),
            this->dimensions()
        )
    );

    Foam::T(result().internalField(), internalField());
    Foam::T(result().boundaryField(), boundaryField());

    return result;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
    >
>
Foam::GeometricField<Type, PatchField, GeoMesh>::component
(
    const direction d
) const
{
    tmp<GeometricField<cmptType, PatchField, GeoMesh> > Component
    (
        new GeometricField<cmptType, PatchField, GeoMesh>
        (
            IOobject
            (
                this->name() + ".component(" + Foam::name(d) + ')',
                this->instance(),
                this->db()
            ),
            this->mesh(),
            this->dimensions()
        )
    );

    Foam::component(Component().internalField(), internalField(), d);
    Foam::component(Component().boundaryField(), boundaryField(), d);

    return Component;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::replace
(
    const direction d,
    const GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
     >& gcf
)
{
    internalField().replace(d, gcf.internalField());
    boundaryField().replace(d, gcf.boundaryField());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::replace
(
    const direction d,
    const dimensioned<cmptType>& ds
)
{
    internalField().replace(d, ds.value());
    boundaryField().replace(d, ds.value());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::max
(
    const dimensioned<Type>& dt
)
{
    Foam::max(internalField(), internalField(), dt.value());
    Foam::max(boundaryField(), boundaryField(), dt.value());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::min
(
    const dimensioned<Type>& dt
)
{
    Foam::min(internalField(), internalField(), dt.value());
    Foam::min(boundaryField(), boundaryField(), dt.value());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::negate()
{
    internalField().negate();
    boundaryField().negate();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator=
(
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
{
    if (this == &gf)
    {
        FatalErrorIn
        (
            "GeometricField<Type, PatchField, GeoMesh>::operator="
            "(const GeometricField<Type, PatchField, GeoMesh>&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    checkField(*this, gf, "=");

    // only equate field contents not ID

    dimensionedInternalField() = gf.dimensionedInternalField();
    boundaryField() = gf.boundaryField();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator=
(
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf
)
{
    if (this == &(tgf()))
    {
        FatalErrorIn
        (
            "GeometricField<Type, PatchField, GeoMesh>::operator="
            "(const tmp<GeometricField<Type, PatchField, GeoMesh> >&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    const GeometricField<Type, PatchField, GeoMesh>& gf = tgf();

    checkField(*this, gf, "=");

    // only equate field contents not ID

    this->dimensions() = gf.dimensions();

    // This is dodgy stuff, don't try it at home.
    internalField().transfer
    (
        const_cast<Field<Type>&>(gf.internalField())
    );

    boundaryField() = gf.boundaryField();

    tgf.clear();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator=
(
    const dimensioned<Type>& dt
)
{
    dimensionedInternalField() = dt;
    boundaryField() = dt.value();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator==
(
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf
)
{
    const GeometricField<Type, PatchField, GeoMesh>& gf = tgf();

    checkField(*this, gf, "==");

    // only equate field contents not ID

    dimensionedInternalField() = gf.dimensionedInternalField();
    boundaryField() == gf.boundaryField();

    tgf.clear();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator==
(
    const dimensioned<Type>& dt
)
{
    dimensionedInternalField() = dt;
    boundaryField() == dt.value();
}


#define COMPUTED_ASSIGNMENT(TYPE, op)                                         \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator op             \
(                                                                             \
    const GeometricField<TYPE, PatchField, GeoMesh>& gf                       \
)                                                                             \
{                                                                             \
    checkField(*this, gf, #op);                                               \
                                                                              \
    dimensionedInternalField() op gf.dimensionedInternalField();              \
    boundaryField() op gf.boundaryField();                                    \
}                                                                             \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator op             \
(                                                                             \
    const tmp<GeometricField<TYPE, PatchField, GeoMesh> >& tgf                \
)                                                                             \
{                                                                             \
    operator op(tgf());                                                       \
    tgf.clear();                                                              \
}                                                                             \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator op             \
(                                                                             \
    const dimensioned<TYPE>& dt                                               \
)                                                                             \
{                                                                             \
    dimensionedInternalField() op dt;                                         \
    boundaryField() op dt.value();                                            \
}

COMPUTED_ASSIGNMENT(Type, +=)
COMPUTED_ASSIGNMENT(Type, -=)
COMPUTED_ASSIGNMENT(scalar, *=)
COMPUTED_ASSIGNMENT(scalar, /=)

#undef COMPUTED_ASSIGNMENT


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
{
    gf.dimensionedInternalField().writeData(os, "internalField");
    os  << nl;
    gf.boundaryField().writeEntry("boundaryField", os);

    // Check state of IOstream
    os.check
    (
        "Ostream& operator<<(Ostream&, "
        "const GeometricField<Type, PatchField, GeoMesh>&)"
    );

    return (os);
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf
)
{
    os << tgf();
    tgf.clear();
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef checkField

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "GeometricBoundaryField.C"
#include "GeometricFieldFunctions.C"

// ************************************************************************* //
