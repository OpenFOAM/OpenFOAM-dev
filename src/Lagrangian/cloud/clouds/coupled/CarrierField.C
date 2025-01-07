/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
#include "coupled.H"

/*---------------------------------------------------------------------------*\
                  Class CloudDerivedField::Functor Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Foam::CarrierField<Type>::Functor
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
class Foam::CarrierField<Type>::Function
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
Foam::CarrierField<Type>::interpolate
(
    const LagrangianModelRef&,
    const LagrangianSubMesh& subMesh
) const
{
    const LagrangianMesh& mesh = subMesh.mesh();

    // Determine whether old-time values are available, and whether or not they
    // should be used
    const bool useOldTime =
        psi().hasStoredOldTimes()
     && psi().nOldTimes(false)
     && mesh.foundObject<LagrangianScalarInternalDynamicField>
        (
            LagrangianMesh::fractionName
        );

    /*
    Info<< "--> Interpolation of carrier field " << name_ << ' '
        << (useOldTime ? "*IS*" : "is *NOT*") << " using old-time values"
        << endl;
    */

    // Construct the interpolation on demand
    if (!psiInterpolationPtr_.valid())
    {
        psiInterpolationPtr_ =
            interpolation<Type>::New
            (
                word(subMesh.mesh().schemes().interpolation(name())),
                psi()
            );
    }

    // Construct the old-time interpolation on demand, if necessary
    if (useOldTime && !psi0InterpolationPtr_.valid())
    {
        psi0InterpolationPtr_ =
            interpolation<Type>::New
            (
                word(subMesh.mesh().schemes().interpolation(name())),
                psi().oldTime()
            );
    }

    // Interpolate the values
    tmp<Field<Type>> tpsic =
        psiInterpolationPtr_->interpolate
        (
            subMesh.sub(mesh.coordinates()),
            subMesh.sub(mesh.celli()),
            subMesh.sub(mesh.facei()),
            subMesh.sub(mesh.faceTrii())
        );

    // If using old-time values, interpolate them also, and combine with the
    // current-time values using the current tracking fraction
    if (useOldTime)
    {
        tmp<Field<Type>> tpsi0c =
            psi0InterpolationPtr_->interpolate
            (
                subMesh.sub(mesh.coordinates()),
                subMesh.sub(mesh.celli()),
                subMesh.sub(mesh.facei()),
                subMesh.sub(mesh.faceTrii())
            );

        SubField<scalar> fraction =
            subMesh.sub
            (
                mesh
               .lookupObject<LagrangianScalarInternalDynamicField>
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
                name() + ':' + Foam::name(subMesh.group()),
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::CarrierField<Type>::CarrierField(const VolField<Type>& psi)
:
    CloudDerivedField<Type>(*this, &CarrierField::interpolate),
    name_(clouds::coupled::carrierName(psi.name())),
    functorPtr_(nullptr),
    tpsi_(psi),
    psiInterpolationPtr_(nullptr),
    psi0InterpolationPtr_(nullptr)
{}


template<class Type>
template<class F>
Foam::CarrierField<Type>::CarrierField(const word& name, const F& f)
:
    CloudDerivedField<Type>(*this, &CarrierField::interpolate),
    name_(name),
    functorPtr_(new Function<F>(name_, f)),
    tpsi_(nullptr),
    psiInterpolationPtr_(nullptr),
    psi0InterpolationPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::VolField<Type>& Foam::CarrierField<Type>::psi() const
{
    if (!tpsi_.valid())
    {
        tpsi_ = functorPtr_->operator()();
    }

    return tpsi_();
}


template<class Type>
const Foam::word& Foam::CarrierField<Type>::name() const
{
    return name_;
}


template<class Type>
void Foam::CarrierField<Type>::reset(const bool predict)
{
    CloudDerivedField<Type>::clear(true);

    if (tpsi_.valid() && tpsi_.isTmp())
    {
        if (predict)
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


// ************************************************************************* //
