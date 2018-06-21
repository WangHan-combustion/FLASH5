# Makefile for source/physics/RadTrans/RadTransTwoMoment

RadTrans += RadTrans_data.o

#  These live in the thornado repo itself

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
RadTrans += TwoMoment_ClosureModule.o
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

RadTrans.o : ProgramHeaderModule.o FluidFieldsModule.o RadiationFieldsModule.o ThornadoInitializationModule.o TimeSteppingModule_Castro.o UnitsModule.o

ProgramHeaderModule.o: \
  KindModule.o

QuadratureModule.o: \
  KindModule.o
ReferenceElementModuleX.o: \
  KindModule.o \
  QuadratureModule.o \
  ProgramHeaderModule.o
ReferenceElementModuleE.o: \
  KindModule.o \
  QuadratureModule.o \
  ProgramHeaderModule.o
ReferenceElementModule.o: \
  KindModule.o \
  QuadratureModule.o \
  ProgramHeaderModule.o
PolynomialBasisModule_Lagrange.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  QuadratureModule.o \
  ReferenceElementModule.o
ReferenceElementModuleX_Lagrange.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  PolynomialBasisModule_Lagrange.o
ReferenceElementModuleE_Lagrange.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleE.o \
  PolynomialBasisModule_Lagrange.o
ReferenceElementModule_Lagrange.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleE.o \
  PolynomialBasisModule_Lagrange.o

MeshModule.o: \
  KindModule.o \
  QuadratureModule.o

GeometryFieldsModule.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  UnitsModule.o
GeometryComputationModule.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  MeshModule.o \
  GeometryFieldsModule.o
GeometryFieldsModuleE.o: \
  KindModule.o \
  ProgramHeaderModule.o
GeometryComputationModuleE.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleE.o \
  ReferenceElementModuleE_Lagrange.o \
  MeshModule.o \
  GeometryFieldsModuleE.o
  
FluidFieldsModule.o: \
  KindModule.o \
  ProgramHeaderModule.o
RadiationFieldsModule.o: \
  KindModule.o \
  ProgramHeaderModule.o

TwoMoment_ClosureModule.o: \
  KindModule.o
TwoMoment_UtilitiesModule.o: \
  KindModule.o
TwoMoment_BoundaryConditionsModule.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModule.o \
  RadiationFieldsModule.o
TwoMoment_PositivityLimiterModule.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModule.o \
  ReferenceElementModule_Lagrange.o \
  RadiationFieldsModule.o
TwoMoment_DiscretizationModule_Streaming.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleX_Lagrange.o \
  ReferenceElementModule.o \
  ReferenceElementModule_Lagrange.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  GeometryComputationModule.o \
  GeometryFieldsModuleE.o \
  RadiationFieldsModule.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_BoundaryConditionsModule.o \
  TwoMoment_UtilitiesModule.o \
  TwoMoment_PositivityLimiterModule.o

PhysicalConstantsModule.o: \
  KindModule.o
UnitsModule.o: \
  KindModule.o \
  PhysicalConstantsModule.o

ThornadoInitializationModule.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  QuadratureModule.o \
  ReferenceElementModuleX.o \
  ReferenceElementModuleE.o \
  ReferenceElementModule.o \
  PolynomialBasisModule_Lagrange.o \
  ReferenceElementModuleX_Lagrange.o \
  ReferenceElementModuleE_Lagrange.o \
  ReferenceElementModule_Lagrange.o \
  MeshModule.o \
  GeometryFieldsModule.o \
  FluidFieldsModule.o \
  GeometryFieldsModuleE.o \
  GeometryComputationModuleE.o \
  RadiationFieldsModule.o \
  TwoMoment_ClosureModule.o \
  TwoMoment_PositivityLimiterModule.o
TimeSteppingModule_Castro.o: \
  KindModule.o \
  ProgramHeaderModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o \
  TwoMoment_DiscretizationModule_Streaming.o \
  DiscretizationModule_Dummy.o
DiscretizationModule_Dummy.o: \
  KindModule.o \
  FluidFieldsModule.o \
  RadiationFieldsModule.o
