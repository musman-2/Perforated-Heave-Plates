  const label& patchIDBW = mesh.boundaryMesh().findPatchID("bottomwall");
  
	if( patchIDBW == -1 )
	{
	  FatalError << "patch not found!" << exit(FatalError);
	}
 
	fixedValueFvPatchField<vector>& VBW =
	    refCast<fixedValueFvPatchField<vector>>(U.boundaryFieldRef()[patchIDBW]);

  scalar tt = runTime.value();
  VBW == vectorField( VBW.size() , vector( 0 , 0 , -Amp*Omega*Foam::sin( Omega*tt ) ) );
  
  forAll( accSource , ii )
  {
  accSource[ii] = vector( 0 , 0 , Amp*Omega*Omega*Foam::cos( Omega*tt )  );
  }