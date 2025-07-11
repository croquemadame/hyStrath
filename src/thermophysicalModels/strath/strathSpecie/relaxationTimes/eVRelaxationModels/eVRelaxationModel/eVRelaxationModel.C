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

#include "eVRelaxationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(eVRelaxationModel, 0);
    defineRunTimeSelectionTable(eVRelaxationModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eVRelaxationModel::eVRelaxationModel
(
    const word& name1,
    const label& lname1,
    const dictionary& dict1,
    const dictionary& dict2,
    const volScalarField& p,
    const volScalarField& Tv
)
:
    dict1_(dict1),
    dict2_(dict2),
    name1_(name1),
    lname1_(lname1),
    p_(p),
    Tv_(Tv),
    eVOverwriteDefault_(readBool(dict1_.subDict("thermalRelaxationModels").subDict("eV").lookup("overwriteDefault"))),
    eVSpeciesDependent_(readBool(dict1_.subDict("thermalRelaxationModels").subDict("eV").lookup("speciesDependent")))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::eVRelaxationModel> Foam::eVRelaxationModel::New
(
    const word& name1,
    const label& lname1,
    const dictionary& dict1,
    const dictionary& dict2,
    const volScalarField& p,
    const volScalarField& Tv
)
{
    word eVRelaxationModelTypeName(dict1.subDict("thermalRelaxationModels").subDict("eV").lookup("model"));

    //dictionaryConstructorTable::iterator cstrIter =
    //    dictionaryConstructorTablePtr_->find(eVRelaxationModelTypeName);
    auto* ctorPtr = dictionaryConstructorTable(eVRelaxationModelTypeName);
    //if (cstrIter == dictionaryConstructorTablePtr_->end())
    if(!ctorPtr)
    {
        FatalErrorIn
        (
            "eVRelaxationModel::New(const volVectorField&, "
            "const surfaceScalarField&)"
        )   << "Unknown eVRelaxationModel type "
            << eVRelaxationModelTypeName << endl << endl
            << "Valid  eVRelaxationModels are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<eVRelaxationModel>
        (ctorPtr(name1, lname1, dict1, dict2, p, Tv));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
