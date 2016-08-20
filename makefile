TARGET = ppp

CODFIL = Main.cpp CFullVec.cpp CHlld.cpp outputPark.cpp setInital.cpp
TRANS = Main.o CFullVec.o CHlld.o outputPark.o setInital.cpp
INCLUD = Param.h CFullVec.h CHlld.h outputPark.h

CC = g++ -fopenmp
KEY = -c


.PHONY: all clear cldir crdir

all: $(TARGET)

clear: 
	rm -rf $(TARGET) *.o *.dat
cldir:
	rm -rf results
crdir:
	mkdir results
	mkdir results/0_Ro results/1_Vx results/2_Vy results/3_Vz 
	mkdir results/4_Bx results/5_By results/6_Bz results/7_E results/8_P results/9_Ptot results/10_div
	#mkdir single_results


Main.o: Main.cpp $(INCLUD)
	$(CC) $(KEY) -o Main.o Main.cpp	

outputPark.o: outputPark.cpp $(INCLUD)
	$(CC) $(KEY) -o outputPark.o outputPark.cpp

CHlld.o: CHlld.cpp $(INCLUD)
	$(CC) $(KEY) -o CHlld.o CHlld.cpp

CFullVec.o: CFullVec.cpp $(INCLUD)
	$(CC) $(KEY) -o CFullVec.o CFullVec.cpp

setInital.o: setInital.cpp $(INCLUD)
	$(CC) $(KEY) -o setInital.o setInital.cpp


$(TARGET): $(TRANS) $(INCLUD)
	$(CC) -o $(TARGET) $(TRANS)