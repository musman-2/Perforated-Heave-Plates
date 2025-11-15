Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField Uss
(
	IOobject
	(
	    "Uss",
	    runTime.timeName(),
	    mesh,
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
	),
	U
); 


volVectorField accSource
(
    IOobject
    (
        "accSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("accSource", dimensionSet(0, 1, -2, 0, 0, 0, 0), vector(0, 0, 0) )
);

#include "createPhi.H"

IOdictionary bcDict
(
    IOobject
    (
        "stateD",    
        runTime.constant(),     
        mesh,                   
        IOobject::MUST_READ_IF_MODIFIED,    
        IOobject::NO_WRITE      
    )
); 

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"

scalar Amp(readScalar(bcDict.lookup("Amplitude")));
scalar Omega(readScalar(bcDict.lookup("Omega")));
