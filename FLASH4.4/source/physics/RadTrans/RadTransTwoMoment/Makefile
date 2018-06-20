# Makefile for source/physics/RadTrans/RadTransTwoMoment

RadTrans += RadTrans_data.o

RadTrans.o : ProgramHeaderModule.o FluidFieldsModule.o RadiationFieldsModule.o ThornadoInitializationModule.o TimeSteppingModule_Castro.o UnitsModule.o

#  These live in the thornado repo itself
  VPATH += ~/thornado/trunk/SandBox/Interface_CASTRO
  VPATH += ~/thornado/trunk/Modules/Numbers
  VPATH += ~/thornado/trunk/Modules/ProgramHeader
  VPATH += ~/thornado/trunk/Modules/Library
  VPATH += ~/thornado/trunk/Modules/Mesh
  VPATH += ~/thornado/trunk/Modules/Fields
  VPATH += ~/thornado/trunk/Modules/Geometry
  VPATH += ~/thornado/trunk/Modules/TwoMoment
  VPATH += ~/thornado/trunk/Modules/Units

# In Modules/Numbers
  RadTrans += KindModule.o

# In Modules/ProgramHeader 
  RadTrans += ProgramHeaderModule.o

# In Modules/Library
  RadTrans += QuadratureModule.o
  RadTrans += ReferenceElementModuleE.o
  RadTrans += ReferenceElementModuleX.o
  RadTrans += ReferenceElementModule.o
  RadTrans += PolynomialBasisModule_Lagrange.o
  RadTrans += ReferenceElementModuleX_Lagrange.o
  RadTrans += ReferenceElementModuleE_Lagrange.o
  RadTrans += ReferenceElementModule_Lagrange.o

# In Modules/Mesh
  RadTrans += MeshModule.o

# In Modules/Geometry
  RadTrans += GeometryFieldsModuleE.o	    
  RadTrans += GeometryFieldsModule.o	   
  RadTrans += GeometryComputationModule.o
  RadTrans += GeometryComputationModuleE.o

# In Modules/Fields
  RadTrans += FluidFieldsModule.o
  RadTrans += RadiationFieldsModule.o

# In Modules/TwoMoment
  ca_F90EXE_sources += TwoMoment_ClosureModule.o
  RadTrans += TwoMoment_PositivityLimiterModule.o
  RadTrans += TwoMoment_DiscretizationModule_Streaming.o
  RadTrans += TwoMoment_BoundaryConditionsModule.o
  RadTrans += TwoMoment_UtilitiesModule.o

# In Modules/Units
  RadTrans += PhysicalConstantsModule.o
  RadTrans += UnitsModule.o

# In SandBox/Interface_CASTRO
  RadTrans += ThornadoInitializationModule.o
  RadTrans += TimeSteppingModule_Castro.o
  RadTrans += DiscretizationModule_Dummy.o
