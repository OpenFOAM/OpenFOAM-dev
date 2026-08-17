/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "fv.H"
#include "objectRegistry.H"
#include "fvMesh.H"
#include "extrapolatedCalculatedFvPatchField.H"
#include "solution.H"

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fv::gradScheme<Type>> Foam::fv::gradScheme<Type>::New
(
    const fvMesh& mesh,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        InfoInFunction << "Constructing gradScheme<Type>" << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Grad scheme not specified" << endl << endl
            << "Valid grad schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(schemeName);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Unknown grad scheme " << schemeName << nl << nl
            << "Valid grad schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::gradScheme<Type>::~gradScheme()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::VolField<typename Foam::outerProduct<Foam::vector, Type>::type>
>
Foam::fv::gradScheme<Type>::fvcGrad_
(
    const VolField<Type>& vf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vf.mesh();

    tmp<VolField<GradType>> tGrad
    (
        VolField<GradType>::New
        (
            name,
            mesh,
            vf.dimensions()/dimensions::length,
            extrapolatedCalculatedFvPatchField<GradType>::typeName
        )
    );
    VolField<GradType>& grad = tGrad.ref();

    calcGrad(grad.internalFieldRef(), vf);

    // Correct the boundary conditions
    grad.correctBoundaryConditions();
    correctBoundaryConditions(vf, grad);

    return tGrad;
}


template<class Type>
void Foam::fv::gradScheme<Type>::correctBoundaryConditions
(
    const VolField<Type>& vf,
    VolField<typename outerProduct<vector, Type>::type>& gGrad
)
{
    typename VolField<typename outerProduct<vector, Type>::type>::Boundary&
        gGradbf = gGrad.boundaryFieldRef();

    forAll(vf.boundaryField(), patchi)
    {
        if (!vf.boundaryField()[patchi].coupled())
        {
            const vectorField n
            (
                vf.mesh().Sf().boundaryField()[patchi]
              / vf.mesh().magSf().boundaryField()[patchi]
            );

            gGradbf[patchi] += n *
            (
                vf.boundaryField()[patchi].snGrad()
              - (n & gGradbf[patchi])
            );
        }
    }
}


template<class Type>
Foam::tmp
<
    Foam::VolInternalField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type
    >
>
Foam::fv::gradScheme<Type>::fviGrad
(
    const VolField<Type>& vf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vf.mesh();

    tmp<VolInternalField<GradType>> tGrad
    (
        VolInternalField<GradType>::New
        (
            name,
            mesh,
            vf.dimensions()/dimensions::length
        )
    );

    calcGrad(tGrad.ref(), vf);

    return tGrad;
}


template<class Type>
Foam::tmp
<
    Foam::VolInternalField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type
    >
>
Foam::fv::gradScheme<Type>::fviGrad
(
    const VolField<Type>& vf
) const
{
    return fviGrad(vf, "grad(" + vf.name() + ')');
}


template<class Type>
Foam::tmp
<
    Foam::VolInternalField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type
    >
>
Foam::fv::gradScheme<Type>::fviGrad
(
    const tmp<VolField<Type>>& tvf
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<VolInternalField<GradType>> tgrad = fviGrad(tvf());
    tvf.clear();
    return tgrad;
}


template<class Type>
Foam::tmp
<
    Foam::VolField<typename Foam::outerProduct<Foam::vector, Type>::type>
>
Foam::fv::gradScheme<Type>::fvcGrad
(
    const VolField<Type>& vf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    if (!this->mesh().changing() && this->mesh().solution().cache(name))
    {
        if
        (
            !mesh().objectRegistry::template
            foundObject<VolField<GradType>>(name)
        )
        {
            solution::cachePrintMessage("Calculating and caching", name, vf);
            tmp<VolField<GradType>> tgGrad = fvcGrad_(vf, name);
            regIOobject::store(tgGrad.ptr());
        }

        solution::cachePrintMessage("Retrieving", name, vf);
        VolField<GradType>& gGrad =
            mesh().objectRegistry::template lookupObjectRef<VolField<GradType>>
            (
                name
            );

        if (gGrad.upToDate(vf))
        {
            return gGrad;
        }
        else
        {
            solution::cachePrintMessage("Deleting", name, vf);
            gGrad.release();
            delete &gGrad;

            solution::cachePrintMessage("Recalculating", name, vf);
            tmp<VolField<GradType>> tgGrad = fvcGrad_(vf, name);

            solution::cachePrintMessage("Storing", name, vf);
            regIOobject::store(tgGrad.ptr());
            VolField<GradType>& gGrad =
                mesh().objectRegistry::template
                lookupObjectRef<VolField<GradType>>
                (
                    name
                );

            return gGrad;
        }
    }
    else
    {
        if
        (
            mesh().objectRegistry::template
            foundObject<VolField<GradType>>(name)
        )
        {
            VolField<GradType>& gGrad =
                mesh().objectRegistry::template
                lookupObjectRef<VolField<GradType>>
                (
                    name
                );

            if (gGrad.ownedByRegistry())
            {
                solution::cachePrintMessage("Deleting", name, vf);
                gGrad.release();
                delete &gGrad;
            }
        }

        solution::cachePrintMessage("Calculating", name, vf);
        return fvcGrad_(vf, name);
    }
}


template<class Type>
Foam::tmp
<
    Foam::VolField<typename Foam::outerProduct<Foam::vector, Type>::type>
>
Foam::fv::gradScheme<Type>::fvcGrad
(
    const VolField<Type>& vf
) const
{
    return fvcGrad(vf, "grad(" + vf.name() + ')');
}


template<class Type>
Foam::tmp
<
    Foam::VolField<typename Foam::outerProduct<Foam::vector, Type>::type>
>
Foam::fv::gradScheme<Type>::fvcGrad
(
    const tmp<VolField<Type>>& tvf
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<VolField<GradType>> tgrad = fvcGrad(tvf());
    tvf.clear();
    return tgrad;
}


// ************************************************************************* //
