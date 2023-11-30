TARGET=macro.exe
all: macro.exe

$(TARGET) : macro.C
		g++ -c `root-config --cflags --libs`  macro.C
		g++ -o $(TARGET)  macro.o myDictionary.o `root-config --cflags --libs`

clean :
	    rm $(TARGET)
