/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

#include "CarrierField.H"
#include "carried.H"
#include "volFields.H"

/*---------------------------------------------------------------------------*\
                  Class CloudDerivedField::Functor Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Foam::CarrierFieldBase<Type>::Functor
{
public:

    // Constructors

        //- Construct null
        Functor()
        {}


    //- Destructor
    virtual ~Functor()
    {}


    // Member Operators

        //- Evaluate the field
        virtual tmp<VolField<Type>> operator()() const = 0;
};


/*---------------------------------------------------------------------------*\
                 Class CloudDerivedField::Function Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
template<class F>
class Foam::CarrierFieldBase<Type>::Function
:
    public Functor
{
    // Private Member Data

        //- The name
        const word& name_;

        //- The function
        F f_;


public:

    // Constructors

        //- Construct from a name and a function
        Function(const word& name, const F& f)
        :
            Functor(),
            name_(name),
            f_(f)
        {}


    //- Destructor
    virtual ~Function()
    {}


    // Member Operators

        //- Evaluate the field
        virtual tmp<VolField<Type>> operator()() const
        {
            tmp<VolField<Type>> tpsi = f_();
            return tpsi.isTmp() ? VolField<Type>::New(name_, tpsi()) : tpsi;
        }
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::CarrierFieldBase<Type>::interpolate
(
    const LagrangianModelRef&,
    const LagrangianSubMesh& subMesh
) const
{
    const LagrangianMesh& mesh = subMesh.mesh();

    tmp<Field<Type>> tpsic =
        psiInterpolation(mesh).interpolate
        (
            subMesh.sub(mesh.coordinates()),
            subMesh.sub(mesh.celli()),
            subMesh.sub(mesh.facei()),
            subMesh.sub(mesh.faceTrii())
        );

    // If using old-time values, interpolate them also, and combine with the
    // current-time values using the current tracking fraction
    if (useOldTime(mesh))
    {
        tmp<Field<Type>> tpsi0c =
            psi0Interpolation(mesh).interpolate
            (
                subMesh.sub(mesh.coordinates()),
                subMesh.sub(mesh.celli()),
                subMesh.sub(mesh.facei()),
                subMesh.sub(mesh.faceTrii())
            );

        const SubField<scalar> fraction =
            subMesh.sub
            (
                mesh.lookupObject<LagrangianInternalScalarDynamicField>
                (
                    LagrangianMesh::fractionName
                ).primitiveField()
            );

        tpsic = (1 - fraction)*tpsi0c + fraction*tpsic;
    }

    // Build the dimensioned field and return
    return tmp<LagrangianSubField<Type>>
    (
        new LagrangianSubField<Type>
        (
            IOobject
            (
                subMesh.sub(this->name_),
                mesh.time().name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            subMesh,
            psi().dimensions(),
            tpsic
        )
    );
}


template<class Type>
Foam::tmp
<
    Foam::LagrangianSubField
    <
        typename Foam::CarrierFieldGradBase<Type>::GradType
    >
>
Foam::CarrierFieldGradBase<Type>::interpolateGrad
(
    const LagrangianModelRef&,
    const LagrangianSubMesh& subMesh
) const
{
    const LagrangianMesh& mesh = subMesh.mesh();

    // Interpolate the values
    tmp<Field<GradType>> tpsic =
        this->psiInterpolation(mesh).interpolateGrad
        (
            subMesh.sub(mesh.coordinates()),
            subMesh.sub(mesh.celli()),
            subMesh.sub(mesh.facei()),
            subMesh.sub(mesh.faceTrii())
        );

    // If using old-time values, interpolate them also, and combine with the
    // current-time values using the current tracking fraction
    if (this->useOldTime(mesh))
    {
        tmp<Field<GradType>> tpsi0c =
            this->psi0Interpolation(mesh).interpolateGrad
            (
                subMesh.sub(mesh.coordinates()),
                subMesh.sub(mesh.celli()),
                subMesh.sub(mesh.facei()),
                subMesh.sub(mesh.faceTrii())
            );

        const SubField<scalar> fraction =
            subMesh.sub
            (
                mesh
               .lookupObject<LagrangianInternalScalarDynamicField>
                (
                    LagrangianMesh::fractionName
                ).primitiveField()
            );

        tpsic = (1 - fraction)*tpsi0c + fraction*tpsic;
    }

    // Build the dimensioned field and return
    return tmp<LagrangianSubField<GradType>>
    (
        new LagrangianSubField<GradType>
        (
            IOobject
            (
                subMesh.sub("grad(" + this->name_ + ')'),
                mesh.time().name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            subMesh,
            this->psi().dimensions()/dimLength,
            tpsic
        )
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::CarrierFieldBase<Type>::useOldTime(const LagrangianMesh& mesh) const
{
    return
        psi().hasStoredOldTimes()
     && psi().nOldTimes(false)
     && mesh.foundObject<LagrangianInternalScalarDynamicField>
        (
            LagrangianMesh::fractionName
        );
}


template<class Type>
const Foam::interpolation<Type>& Foam::CarrierFieldBase<Type>::psiInterpolation
(
    const LagrangianMesh& mesh
) const
{
    if (!psiInterpolationPtr_.valid())
    {
        psiInterpolationPtr_ =
            interpolation<Type>::New
            (
                word(mesh.schemes().interpolation(this->name_)),
                psi()
            );
    }

    return psiInterpolationPtr_();
}


template<class Type>
const Foam::interpolation<Type>& Foam::CarrierFieldBase<Type>::psi0Interpolation
(
    const LagrangianMesh& mesh
) const
{
    if (!psi0InterpolationPtr_.valid())
    {
        psi0InterpolationPtr_ =
            interpolation<Type>::New
            (
                word(mesh.schemes().interpolation(this->name_)),
                psi().oldTime()
            );
    }

    return psi0InterpolationPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::CarrierFieldBase<Type>::CarrierFieldBase(const VolField<Type>& psi)
:
    CloudDerivedField<Type>
    (
        clouds::carried::nameToCarrierName(psi.name()),
        *this,
        &CarrierFieldBase::interpolate
    ),
    functorPtr_(nullptr),
    tpsi_(psi),
    psiInterpolationPtr_(nullptr),
    psi0InterpolationPtr_(nullptr)
{}


template<class Type>
Foam::CarrierFieldBase<Type>::CarrierFieldBase
(
    const word& name,
    const VolField<Type>& psi
)
:
    CloudDerivedField<Type>(name, *this, &CarrierFieldBase::interpolate),
    functorPtr_(nullptr),
    tpsi_(psi),
    psiInterpolationPtr_(nullptr),
    psi0InterpolationPtr_(nullptr)
{}


template<class Type>
template<class F>
Foam::CarrierFieldBase<Type>::CarrierFieldBase(const word& name, const F& f)
:
    CloudDerivedField<Type>(name, *this, &CarrierFieldBase::interpolate),
    functorPtr_(new Function<F>(this->name_, f)),
    tpsi_(nullptr),
    psiInterpolationPtr_(nullptr),
    psi0InterpolationPtr_(nullptr)
{}


template<class Type>
template<class ... Args>
Foam::CarrierFieldGradBase<Type>::CarrierFieldGradBase(const Args& ... args)
:
    CarrierFieldBase<Type>(args ...),
    grad(*this, &CarrierFieldGradBase::interpolateGrad)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::VolField<Type>& Foam::CarrierFieldBase<Type>::psi() const
{
    if (!tpsi_.valid())
    {
        tpsi_ = functorPtr_->operator()();
    }

    return tpsi_();
}


template<class Type>
void Foam::CarrierFieldBase<Type>::reset(const bool initial)
{
    CloudDerivedField<Type>::clear(true);

    if (tpsi_.valid() && tpsi_.isTmp())
    {
        if (initial)
        {
            tpsi_.ref().nullOldestTime();
        }
        else
        {
            tpsi_.ref().oldTime();
            tpsi_.ref() = functorPtr_->operator()();
        }
    }

    psiInterpolationPtr_.clear();
    psi0InterpolationPtr_.clear();
}


template<class Type>
void Foam::CarrierFieldGradBase<Type>::reset(const bool initial)
{
    CarrierFieldBase<Type>::reset(initial);

    grad.clear(true);
}


// ************************************************************************* //
