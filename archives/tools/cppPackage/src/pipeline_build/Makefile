SHELL = /bin/sh

################################################################################
# Principal Investigator : Joel Graber
# Author : Lucie Ndzana Hutchins
# Title  : Scientific Software Engineer
# Company : Computational Science and System Biology, The Jackson Laboratory 
# 
# The main makefile that generates the library (libcpp_core.a) needed to build
# the nnmf, and the executables of all the nnmf's supporting programs
#  -D_FILE_OFFSET_BITS=64
#################################################################################
#patterns
#%.o: %.cpp %.h
#        gcc $(PRODUCTFLAGS) -c $<
#%.do: %.cpp %.h
#        gcc $(DEBUGFLAGS) -c $<
#%.io: %.cpp %.h
#        icpc $(IPRODUCTFLAGS) -DINTEL -o $@ -c $<
#%.h:

CC= g++

MAINLIB = PosWord.o arg.o basicSeq.o \
	fastaSeq.o filter.o logFile.o  \
	util.o vectorF.o wordCount.o wordStat.o 

OBJS = PosWord.o align.o arg.o basicSeq.o bayesWord.o direct.o \
	fastaSeq.o filter.o genome.o logFile.o matrixF.o \
	util.o vectorF.o wCount.o wordCount.o wordStat.o 
	
INC = -I/usr/include
INC1 = -I/usr/local/include
INC2 = -I ./lib
CFLAG = -c -w -O -Wno-deprecated  
CCFLAGS =$(CFLAG) $(INC) $(INC1) $(INC2)

LFLAG4 = -L/ubuntu/usr/lib  -lgsl -lgslcblas
LFLAG1 =  -L/sw/lib -lgsl -lgslcblas 
LFLAG2 =  -L/usr/local/lib
LFLAG3 =  -L/usr/lib
LIBCORE = ./lib

all: libcpp_core.a fa2Stats WindowCount  nmfSmoothAndPseudoCounts nmfSortMatrix nmfBuildMotifs nmfMotifsToModels pwmToExamples nmfWebPage 

libcpp_core.a: $(OBJS)
	ar rc libcpp_core.a *.o;

	
fa2Stats: ${MAINLIB}  fa2Stats.o
	${CC}   -o fa2Stats ${MAINLIB} fa2Stats.o ${LFLAG1}  ${LFLAG2}  ${LFLAG3} ${LFLAG4}
	
WindowCount: ${MAINLIB}  WindowCount.o
	${CC}   -o WindowCount ${MAINLIB} WindowCount.o ${LFLAG1}  ${LFLAG2}  ${LFLAG3} ${LFLAG4} 
	
nmfSmoothAndPseudoCounts: ${MAINLIB} nmfSmoothAndPseudoCounts.o
	${CC}   -o nmfSmoothAndPseudoCounts  ${MAINLIB} nmfSmoothAndPseudoCounts.o ${LFLAG1}  ${LFLAG2}  ${LFLAG3} ${LFLAG4}

nmfSortMatrix: ${MAINLIB} nmfSortMatrix.o
	${CC}    -o  nmfSortMatrix  ${MAINLIB} nmfSortMatrix.o ${LFLAG1}  ${LFLAG2}  ${LFLAG3} ${LFLAG4}

nmfBuildMotifs: ${MAINLIB}  buildMotifs.o
	${CC}   -o nmfBuildMotifs ${MAINLIB} buildMotifs.o ${LFLAG1}  ${LFLAG2}  ${LFLAG3} ${LFLAG4}

nmfMotifsToModels: ${MAINLIB} nmfMotifsToModels.o
	${CC}   -o nmfMotifsToModels  ${MAINLIB} nmfMotifsToModels.o ${LFLAG1}  ${LFLAG2}  ${LFLAG3} ${LFLAG4}

pwmToExamples: ${MAINLIB} pwmToExamples.o
	${CC}  -o pwmToExamples  ${MAINLIB} pwmToExamples.o ${LFLAG1}  ${LFLAG2}  ${LFLAG3} ${LFLAG4}

nmfWebPage: ${MAINLIB} nmfWebPage.o
	${CC}  -o nmfWebPage  ${MAINLIB} nmfWebPage.o ${LFLAG1}  ${LFLAG2}  ${LFLAG3} ${LFLAG4}


fa2Stats.o: ${LIBCORE}/fa2Stats.cpp ${LIBCORE}/basicSeq.h ${LIBCORE}/arg.h ${LIBCORE}/wordCount.h ${LIBCORE}/fastaSeq.h ${LIBCORE}/logFile.h
	${CC} ${CCFLAGS} ${LIBCORE}/fa2Stats.cpp

WindowCount.o: ${LIBCORE}/WindowCount.cpp ${LIBCORE}/basicSeq.h ${LIBCORE}/arg.h ${LIBCORE}/PosWord.h ${LIBCORE}/logFile.h
	${CC} ${CCFLAGS} ${LIBCORE}/WindowCount.cpp

nmfSmoothAndPseudoCounts.o : ${LIBCORE}/nmfSmoothAndPseudoCounts.cpp  ${LIBCORE}/arg.h ${LIBCORE}/util.h ${LIBCORE}/logFile.h
	${CC} ${CCFLAGS} ${LIBCORE}/nmfSmoothAndPseudoCounts.cpp
	
nmfSortMatrix.o : ${LIBCORE}/nmfSortMatrix.cpp  ${LIBCORE}/arg.h ${LIBCORE}/util.h ${LIBCORE}/logFile.h
	${CC} ${CCFLAGS} ${LIBCORE}/nmfSortMatrix.cpp
	
logFile.o: ${LIBCORE}/logFile.cpp ${LIBCORE}/logFile.h
	${CC} ${CCFLAGS} ${LIBCORE}/logFile.cpp

arg.o:	${LIBCORE}/arg.h ${LIBCORE}/arg.cpp ${LIBCORE}/logFile.h
	${CC} ${CCFLAGS} ${LIBCORE}/arg.cpp
	
basicSeq.o: ${LIBCORE}/basicSeq.h ${LIBCORE}/basicSeq.cpp ${LIBCORE}/logFile.h ${LIBCORE}/util.h
	${CC} ${CCFLAGS} ${LIBCORE}/basicSeq.cpp
	
PosWord.o: ${LIBCORE}/PosWord.cpp ${LIBCORE}/PosWord.h ${LIBCORE}/basicSeq.h ${LIBCORE}/logFile.h
	${CC} ${CCFLAGS} ${LIBCORE}/PosWord.cpp

fastaSeq.o: ${LIBCORE}/fastaSeq.cpp ${LIBCORE}/fastaSeq.h
	${CC} ${CCFLAGS} ${LIBCORE}/fastaSeq.cpp	

util.o:	${LIBCORE}/util.h ${LIBCORE}/util.cpp ${LIBCORE}/logFile.h
	${CC} ${CCFLAGS} ${LIBCORE}/util.cpp
	
filter.o: ${LIBCORE}/filter.h ${LIBCORE}/filter.cpp ${LIBCORE}/logFile.h
	${CC} ${CCFLAGS} ${LIBCORE}/filter.cpp

wordCount.o: ${LIBCORE}/wordCount.cpp ${LIBCORE}/wordCount.h ${LIBCORE}/wordStat.h ${LIBCORE}/logFile.h ${LIBCORE}/util.h ${LIBCORE}/basicSeq.h ${LIBCORE}/vectorF.h
	${CC} ${CCFLAGS} ${LIBCORE}/wordCount.cpp

wordStat.o: ${LIBCORE}/wordStat.cpp ${LIBCORE}/wordStat.h ${LIBCORE}/logFile.h ${LIBCORE}/util.h
	${CC} ${CCFLAGS} ${LIBCORE}/wordStat.cpp

buildMotifs.o : ${LIBCORE}/buildMotifs.cpp  ${LIBCORE}/arg.h ${LIBCORE}/util.h ${LIBCORE}/logFile.h ${LIBCORE}/basicSeq.h ${LIBCORE}/vectorF.h
	${CC} ${CCFLAGS} ${LIBCORE}/buildMotifs.cpp
	
nmfMotifsToModels.o : ${LIBCORE}/nmfMotifsToModels.cpp  ${LIBCORE}/arg.h ${LIBCORE}/util.h ${LIBCORE}/logFile.h
	${CC} ${CCFLAGS} ${LIBCORE}/nmfMotifsToModels.cpp
	
pwmToExamples.o : ${LIBCORE}/pwmToExamples.cpp  ${LIBCORE}/arg.h ${LIBCORE}/util.h ${LIBCORE}/logFile.h ${LIBCORE}/basicSeq.h ${LIBCORE}/fastaSeq.h
	${CC} ${CCFLAGS} ${LIBCORE}/pwmToExamples.cpp
	
nmfWebPage.o : ${LIBCORE}/nmfWebPage.cpp  ${LIBCORE}/arg.h ${LIBCORE}/util.h ${LIBCORE}/logFile.h
	${CC} ${CCFLAGS} ${LIBCORE}/nmfWebPage.cpp

vectorF.o: ${LIBCORE}/vectorF.cpp ${LIBCORE}/vectorF.h  ${LIBCORE}/logFile.h ${LIBCORE}/util.h
	${CC} ${CCFLAGS} ${LIBCORE}/vectorF.cpp

genome.o: ${LIBCORE}/genome.h ${LIBCORE}/genome.cpp ${LIBCORE}/logFile.h ${LIBCORE}/basicSeq.h
	${CC} ${CCFLAGS} ${LIBCORE}/genome.cpp

matrixF.o:	${LIBCORE}/matrixF.cpp ${LIBCORE}/matrixF.h ${LIBCORE}/logFile.h
	${CC} ${CCFLAGS} ${LIBCORE}/matrixF.cpp	

align.o: ${LIBCORE}/align.cpp ${LIBCORE}/align.h ${LIBCORE}/logFile.h ${LIBCORE}/basicSeq.h
	${CC} ${CCFLAGS} ${LIBCORE}/align.cpp

bayesWord.o: ${LIBCORE}/bayesWord.cpp ${LIBCORE}/bayesWord.h ${LIBCORE}/basicSeq.h ${LIBCORE}/logFile.h
	${CC} ${CCFLAGS} ${LIBCORE}/bayesWord.cpp	

direct.o: ${LIBCORE}/direct.cpp ${LIBCORE}/direct.h ${LIBCORE}/basicSeq.h ${LIBCORE}/logFile.h
	${CC} ${CCFLAGS} ${LIBCORE}/direct.cpp	

wCount.o: ${LIBCORE}/wCount.cpp ${LIBCORE}/wCount.h
	${CC} ${CCFLAGS} ${LIBCORE}/wCount.cpp
		
clean: rm  *.o
