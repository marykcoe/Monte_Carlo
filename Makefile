CC = g++

OBJECTS =  $(OBJDIR)particle.o $(OBJDIR)cell.o $(OBJDIR)simbox.o $(OBJDIR)run.o \
			$(OBJDIR)main.o $(OBJDIR)energetics.o $(OBJDIR)measures.o $(OBJDIR)insert.o $(OBJDIR)delete.o \
			$(OBJDIR)output.o $(OBJDIR)solute.o $(OBJDIR)input.o $(OBJDIR)test.o

CFLAGS = -std=c++17 -O3 -g

TARGET = MC

OBJDIR = ./Program/obj/
SRCDIR = ./Program/src/

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJECTS)

$(OBJECTS): | $(OBJDIR)

$(OBJDIR)particle.o: $(SRCDIR)particle.cpp $(SRCDIR)particle.h
		$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)cell.o: $(SRCDIR)cell.cpp $(SRCDIR)cell.h
		$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)solute.o: $(SRCDIR)solute.cpp $(SRCDIR)solute.h $(SRCDIR)simbox.h $(SRCDIR)constants.h
		$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)simbox.o: $(SRCDIR)simbox.cpp $(SRCDIR)simbox.h $(SRCDIR)cell.h \
					$(SRCDIR)particle.h $(SRCDIR)constants.h $(SRCDIR)solute.h
		$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)test.o: $(SRCDIR)test.cpp $(SRCDIR)test.h $(SRCDIR)simbox.h $(SRCDIR)run.h 
		$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)input.o: $(SRCDIR)input.cpp $(SRCDIR)input.h $(SRCDIR)simbox.h  $(SRCDIR)test.h
		$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)delete.o: $(SRCDIR)delete.cpp $(SRCDIR)delete.h $(SRCDIR)cell.h $(SRCDIR)particle.h \
					$(SRCDIR)simbox.h $(SRCDIR)energetics.h
		$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)insert.o: $(SRCDIR)insert.cpp $(SRCDIR)insert.h $(SRCDIR)cell.h $(SRCDIR)particle.h \
					$(SRCDIR)simbox.h $(SRCDIR)energetics.h
		$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)measures.o: $(SRCDIR)measures.cpp $(SRCDIR)measures.h $(SRCDIR)cell.h $(SRCDIR)constants.h \
					$(SRCDIR)particle.h $(SRCDIR)simbox.h
		$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)output.o: $(SRCDIR)output.cpp $(SRCDIR)output.h $(SRCDIR)cell.h $(SRCDIR)constants.h \
					$(SRCDIR)particle.h $(SRCDIR)simbox.h
		$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)run.o: $(SRCDIR)run.cpp $(SRCDIR)run.h $(SRCDIR)cell.h $(SRCDIR)constants.h $(SRCDIR)particle.h \
				$(SRCDIR)measures.h $(SRCDIR)simbox.h $(SRCDIR)delete.h $(SRCDIR)insert.h
		$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)energetics.o: $(SRCDIR)energetics.cpp $(SRCDIR)energetics.h $(SRCDIR)cell.h $(SRCDIR)constants.h \
						$(SRCDIR)particle.h $(SRCDIR)simbox.h 
		$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)main.o: $(SRCDIR)main.cpp $(SRCDIR)cell.h $(SRCDIR)constants.h $(SRCDIR)particle.h $(SRCDIR)simbox.h \
				$(SRCDIR)run.h $(SRCDIR)energetics.h $(SRCDIR)output.h $(SRCDIR)test.h
		$(CC) $(CFLAGS) -c $< -o $@  

$(OBJDIR): 
		mkdir -p $(OBJDIR)

clean:
	rm $(TARGET) $(OBJDIR)*.o 
