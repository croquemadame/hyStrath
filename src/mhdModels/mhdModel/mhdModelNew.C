/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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

#include "mhdModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mhd::mhdModel>
Foam::mhd::mhdModel::New
(
    const rho2ReactionThermo& thermo
)
{
    IOobject mhdIO
    (
        "mhdProperties",
        thermo.T().time().constant(),
        thermo.T().mesh(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    );

    word modelType("noMHD");
    bool active = false;
    
    if (mhdIO.typeHeaderOk<IOobject>())
    {
        if (IOdictionary(mhdIO).found("active"))
        {
            IOdictionary(mhdIO).lookup("active") >> active;
            
            if (active && IOdictionary(mhdIO).found("mhdModel"))
            {
                IOdictionary(mhdIO).lookup("mhdModel") >> modelType;
            }
        }
    }
    else
    {
        Info<< "MHD model not active: mhdProperties dict not found" << endl;
    }

    Info<< "Selecting mhdModel " << modelType << nl << endl;

    auto cstrIter =
        thermoConstructorTablePtr_->find(modelType);

    if (cstrIter == thermoConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "mhdModel::New(const volScalarField&)"
        )   << "Unknown mhdModel type "
            << modelType << nl << nl
            << "Valid mhdModel types are:" << nl
            << thermoConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<mhdModel>(cstrIter()(thermo));
}


Foam::autoPtr<Foam::mhd::mhdModel>
Foam::mhd::mhdModel::New
(
    const dictionary& dict,
    const rho2ReactionThermo& thermo
)
{
    bool active(dict.lookupOrDefault<bool>("active", false));
    word modelType("noMHD");
    
    if (active && dict.found("mhdModel"))
    {
        modelType = word(dict.lookup("mhdModel"));
    }

    Info<< "Selecting mhdModel " << modelType << nl << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "mhdModel::New(const dictionary&, const volScalarField&)"
        )   << "Unknown mhdModel type "
            << modelType << nl << nl
            << "Valid mhdModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<mhdModel>(cstrIter()(dict, thermo));
}


// ************************************************************************* //
