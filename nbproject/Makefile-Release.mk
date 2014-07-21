#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=cc
CCC=CC
CXX=CC
FC=f95
AS=as

# Macros
CND_PLATFORM=OracleSolarisStudio-Solaris-x86
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/Panel.o \
	${OBJECTDIR}/src/PotentialFlowElements.o \
	${OBJECTDIR}/src/AerodynamicBody2D.o \
	${OBJECTDIR}/src/testVector.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ndeaero

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ndeaero: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ndeaero ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/src/Panel.o: src/Panel.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 -o ${OBJECTDIR}/src/Panel.o src/Panel.cpp

${OBJECTDIR}/src/PotentialFlowElements.o: src/PotentialFlowElements.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 -o ${OBJECTDIR}/src/PotentialFlowElements.o src/PotentialFlowElements.cpp

${OBJECTDIR}/src/AerodynamicBody2D.o: src/AerodynamicBody2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 -o ${OBJECTDIR}/src/AerodynamicBody2D.o src/AerodynamicBody2D.cpp

${OBJECTDIR}/src/testVector.o: src/testVector.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.cc) -fast -g0 -o ${OBJECTDIR}/src/testVector.o src/testVector.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ndeaero
	${CCADMIN} -clean

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
