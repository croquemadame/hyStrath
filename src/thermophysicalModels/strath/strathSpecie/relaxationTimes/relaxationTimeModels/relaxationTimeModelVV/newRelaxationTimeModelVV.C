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

#include "relaxationTimeModelVV.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::relaxationTimeModelVV>
Foam::relaxationTimeModelVV::New
(
    rho2ReactionThermo& thermo,
    const compressibleTurbulenceModel& turbulence
)
{
    word VVModelName;

    // Enclose the creation of the dictionary to ensure it is deleted
    // before the relaxationTimeModelVV is created otherwise the
    // dictionary is entered in the database twice
    {
        const dictionary& thermo2TModel =
        (
            IFstream
            (
                fileName(thermo.lookup("twoTemperatureDictFile")).expand()
            )()
        );

        word partialVVname =
            word
            (
                thermo2TModel.subDict("thermalRelaxationModels").subDict("VV")
                    .lookup("relaxationType")
            );

        VVModelName = partialVVname +'<' + thermo.partialThermoName() + '>';

        Info<< "Loading the V-V relaxation time model:" << tab
            << partialVVname << "\n" << endl;
    }

    //fvMeshConstructorTable::iterator cstrIter =
    //    fvMeshConstructorTablePtr_->find(VVModelName);
    auto* ctorPtr = fvMeshConstructorTable(VVModelName);


    //if (cstrIter == fvMeshConstructorTablePtr_->end())
    if(!ctorPtr)
    {
        FatalErrorIn
        (
            "VVModel::New(const volVectorField&, "
            "const surfaceScalarField&)"
        )   << "Unknown VVModel type "
            << VVModelName << endl << endl
            << "Valid  VVModels are : " << endl
            << fvMeshConstructorTablePtr_->toc()
            << exit(FatalError);
    }


    return autoPtr<relaxationTimeModelVV>
        (ctorPtr(thermo, turbulence));
}


// ************************************************************************* //
