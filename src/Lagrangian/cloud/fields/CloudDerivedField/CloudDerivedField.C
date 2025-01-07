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

#include "CloudDerivedField.H"
#include "calculatedLagrangianPatchField.H"
#include "Time.H"

/*---------------------------------------------------------------------------*\
                  Class CloudDerivedField::Functor Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Foam::CloudDerivedField<Type>::Functor
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

        //- Evaluate
        virtual tmp<LagrangianSubField<Type>> operator()
        (
            const LagrangianModelRef& model,
            const LagrangianSubMesh& subMesh
        ) const = 0;
};


/*---------------------------------------------------------------------------*\
                 Class CloudDerivedField::Function Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
template<class F>
class Foam::CloudDerivedField<Type>::Function
:
    public Functor
{
    // Private Member Data

        //- The function
        F f_;


public:

    // Constructors

        //- Construct from a function
        Function(const F& f)
        :
            Functor(),
            f_(f)
        {}


    //- Destructor
    virtual ~Function()
    {}


    // Member Operators

        //- Evaluate the field
        virtual tmp<LagrangianSubField<Type>> operator()
        (
            const LagrangianModelRef& model,
            const LagrangianSubMesh& subMesh
        ) const
        {
            return f_(model, subMesh);
        }
};


/*---------------------------------------------------------------------------*\
                  Class CloudDerivedField::Method Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
template<class C>
class Foam::CloudDerivedField<Type>::Method
:
    public Functor
{
    // Private Member Data

        //- The class
        const C& c_;

        //- The method
        tmp<LagrangianSubField<Type>> (C::*m_)
        (
            const LagrangianModelRef&,
            const LagrangianSubMesh&
        ) const;


public:

    // Constructors

        //- Construct from a class and a method
        Method
        (
            const C& c,
            tmp<LagrangianSubField<Type>> (C::*m)
            (
                const LagrangianModelRef&,
                const LagrangianSubMesh&
            ) const
        )
        :
            Functor(),
            c_(c),
            m_(m)
        {}


    //- Destructor
    virtual ~Method()
    {}


    // Member Operators

        //- Evaluate the field
        virtual tmp<LagrangianSubField<Type>> operator()
        (
            const LagrangianModelRef& model,
            const LagrangianSubMesh& subMesh
        ) const
        {
            return (c_.*m_)(model, subMesh);
        }
};


/*---------------------------------------------------------------------------*\
              Class CloudDerivedField::AllFieldToField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Foam::CloudDerivedField<Type>::AllFieldToField
:
    public LagrangianInternalField<Type>
{
public:

    // Constructors

        //- Construct from a sub-all-field reference
        AllFieldToField(const LagrangianSubField<Type>& allField)
        :
            LagrangianInternalField<Type>
            (
                IOobject
                (
                    allField.name(),
                    allField.mesh().mesh().time().name(),
                    allField.mesh().mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                allField.mesh().mesh(),
                allField.dimensions(),
                NullObjectRef<Field<Type>>()
            )
        {
            this->UList<Type>::shallowCopy(allField);
        }


    //- Destructor
    virtual ~AllFieldToField()
    {
        this->UList<Type>::shallowCopy(UList<Type>(nullptr, 0));
    }
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
template<class F>
Foam::CloudDerivedField<Type>::CloudDerivedField(const word& name, const F& f)
:
    name_(name),
    functorPtr_(new Function<F>(f))
{}


template<class Type>
template<class F>
Foam::CloudDerivedField<Type>::CloudDerivedField(const F& f)
:
    name_(NullObjectRef<word>()),
    functorPtr_(new Function<F>(f))
{}


template<class Type>
template<class C>
Foam::CloudDerivedField<Type>::CloudDerivedField
(
    const word& name,
    const C& c,
    tmp<LagrangianSubField<Type>> (C::*m)
    (
        const LagrangianModelRef&,
        const LagrangianSubMesh&
    ) const
)
:
    name_(name),
    functorPtr_(new Method<C>(c, m))
{}


template<class Type>
template<class C>
Foam::CloudDerivedField<Type>::CloudDerivedField
(
    const C& c,
    tmp<LagrangianSubField<Type>> (C::*m)
    (
        const LagrangianModelRef&,
        const LagrangianSubMesh&
    ) const
)
:
    name_(NullObjectRef<word>()),
    functorPtr_(new Method<C>(c, m))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::LagrangianInternalField<Type>>
Foam::CloudDerivedField<Type>::field(const LagrangianMesh& mesh) const
{
    const LagrangianSubField<Type>& allField =
        operator()(LagrangianModelRef(), mesh.subAll());

    tmp<LagrangianInternalField<Type>> tResult =
        LagrangianInternalField<Type>::New
        (
            allField.name(),
            mesh,
            allField.dimensions()
        );

    tResult.ref().primitiveFieldRef() = allField.primitiveField();

    return tResult;
}


template<class Type>
Foam::LagrangianSubSubField<Type>& Foam::CloudDerivedField<Type>::ref
(
    const LagrangianSubMesh& subMesh
) const
{
    // Evaluate and store if it doesn't already exist for the sub-mesh
    if (!psiSubSubPtr_.valid() || psiSubSubMeshIndex_ != subMesh.index())
    {
        psiSubSubPtr_.reset
        (
            new LagrangianSubSubField<Type>(operator()(subMesh))
        );

        psiSubSubUpToDate_ = true;
        psiSubSubMeshIndex_ = subMesh.index();
    }

    // Update the field in-place if the up-to-date flag is not set
    if (!psiSubSubUpToDate_)
    {
        psiSubSubPtr_->UList<Type>::shallowCopy(operator()(subMesh));

        psiSubSubUpToDate_ = true;
    }

    return psiSubSubPtr_();
}


template<class Type>
void Foam::CloudDerivedField<Type>::clear(const bool final)
{
    psiAllPtr_.clear();

    // If this is not the final iteration, then retain the fields to be
    // modified in-place. This is more efficient, and it ensures that old-time
    // fields are maintained.
    if (!final)
    {
        psiSubUpToDate_ = false;
        psiSubSubUpToDate_ = false;
    }
    else
    {
        psiSubPtr_.clear();
        psiSubSubPtr_.clear();
    }
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::LagrangianInternalField<Type>>
Foam::CloudDerivedField<Type>::operator()(const LagrangianMesh& mesh) const
{
    return
        tmp<LagrangianInternalField<Type>>
        (
            new AllFieldToField
            (
                operator()(LagrangianModelRef(), mesh.subAll())
            ),
            true
        );
}


template<class Type>
const Foam::LagrangianSubField<Type>&
Foam::CloudDerivedField<Type>::operator()
(
    const LagrangianSubMesh& subMesh
) const
{
    return operator()(LagrangianModelRef(), subMesh);
}


template<class Type>
const Foam::LagrangianSubField<Type>&
Foam::CloudDerivedField<Type>::operator()
(
    const LagrangianModelRef& model,
    const LagrangianSubMesh& subMesh
) const
{
    // Evaluate and store if it doesn't already exist for the all-mesh
    if (&subMesh == &subMesh.mesh().subAll())
    {
        if (!psiAllPtr_.valid())
        {
            if (notNull(name_))
            {
                psiAllPtr_.reset
                (
                    new Foam::LagrangianSubField<Type>
                    (
                        name_,
                        functorPtr_()(model, subMesh)
                    )
                );
            }
            else
            {
                psiAllPtr_.reset(functorPtr_()(model, subMesh).ptr());
            }
        }

        return psiAllPtr_();
    }

    // Evaluate and store if it doesn't already exist for the sub-mesh
    if (!psiSubPtr_.valid() || psiSubMeshIndex_ != subMesh.index())
    {
        if (notNull(name_))
        {
            psiSubPtr_.reset
            (
                new Foam::LagrangianSubField<Type>
                (
                    name_ + ':' + name(subMesh.group()),
                    functorPtr_()(model, subMesh)
                )
            );
        }
        else
        {
            psiSubPtr_.reset(functorPtr_()(model, subMesh).ptr());
        }

        psiSubUpToDate_ = true;
        psiSubMeshIndex_ = subMesh.index();
    }

    // Update the field in-place if the up-to-date flag is not set
    if (!psiSubUpToDate_)
    {
        psiSubPtr_() = functorPtr_()(model, subMesh);

        psiSubUpToDate_ = true;
    }

    return psiSubPtr_();
}


// ************************************************************************* //
