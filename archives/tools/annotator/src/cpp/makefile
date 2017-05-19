#makefile for gutils project
GCC = g++
#Intel compiler
#GCC = icpc


CFLAG =  -c -w -Wno-deprecated -D_FILE_OFFSET_BITS=64 
# -std=gnu++0x

LDFLAG1 = -L/usr/local/lib

OBJS = genomeTest.o ggenome.o

OBJSS = genomeAnnotator.o ggenome.o

#OBJ = getTxFastaSeqOpt.o ggenome.o

OBJJ = getTxFastaSeq.o ggenome.o

OBJG = getConnectedComponents.o 
#Create a library
#libcpp_core.a: $(OBJS)
#	ar rc libcpp_core.a *.o;
################## object file creation
#OBJS = snpvalidator.o gfile.o gsequence.o
#
#OBJS = getCDS.o ensembldb.o gsequence.o
#
#OBJS =  getFastaSeq.o gsequence.o gfile.o 
#
#OBJS = snpAnnotator.o aminoAcid.o ensembldb.o gfile.o gsequence.o
################### executable files generation
#getCDS :	$(OBJS)
#	$(GCC) -o getCDS $(OBJS) $(LDFLAG1)
#
#getTxFastaSeq:	$(OBJJ)
#	$(GCC) -o getTxFastaSeq $(OBJJ) $(LDFLAG1)

#getTxFastaSeqOpt:	$(OBJ)
#	$(GCC) -o getTxFastaSeqOpt $(OBJ) $(LDFLAG1)
#
#getConnectedComponents:	$(OBJG)
#	$(GCC) -o getConnectedComponents $(OBJG) $(LDFLAG1)

#snpValidator:	$(OBJS)
#	$(GCC) -o snpValidator $(OBJS) $(LDFLAG1)
#
#snpAnnotator:	$(OBJS)
#	$(GCC) -o snpAnnotator $(OBJS) $(LDFLAG1)
###################################################

genomeAnnotator:	$(OBJSS)
	$(GCC) -o genomeAnnotator $(OBJSS) $(LDFLAG1)

###################### dependencies

#getTxFastaSeq.o: getTxFastaSeq.cc
#	$(GCC) -c getTxFastaSeq.cc

#getFastaSeqOpt.o: getTxFastaSeqOpt.cc
#	$(GCC) $(CFLAG) getTxFastaSeqOpt.cc

#getConnectedComponents.o: getConnectedComponents.cc
#	$(GCC) $(CFLAG) getConnectedComponents.cc

#snpvalidator.o: snpvalidator.cc
#	$(GCC) $(CFLAG) snpvalidator.cc
#
#snpAnnotator.o: snpAnnotator.cc
#	$(GCC) $(CFLAG) snpAnnotator.cc

genomeAnnotator.o: genomeAnnotator.cc
	$(GCC) -c genomeAnnotator.cc

ggenome.o:	 ggenome.h  ggenome.cc
	$(GCC) -c ggenome.cc
	

clean:
	rm $(OBJS) 
