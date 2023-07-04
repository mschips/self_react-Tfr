# Load default variables
include ../defaults.mk

### TARGET SPECIFIC RULES ###

### COMPILE COMMAND ###
COMPILE = $(CXX) $(CPPFLAGS) $(INCLUDES) $(CXXFLAGS)
DCOMPILE = $(CXX) $(CPPFLAGS) $(INCLUDES) $(CXXFLAGS_D)

.PHONY: all dirs clean
all:
	$(CMDPRE)$(MAKE) $(BINPATH)/$(BINARY) --no-print-directory

dirs :
	$(MKDIR) $(BUILDPATH) $(BINPATH)

clean:
	$(RM) $(OBJECTS) $(BINPATH)/$(BINARY) $(DEPS) $(OBJECTS_D) $(DEPS_D)

debug:
	$(CMDPRE)$(MAKE) $(BINPATH)/$(BINARY)_D

### OBJECTS ###
$(OBJECTS): $(BUILDPATH)/%.o: $(SRCPATH)/%.cpp | dirs
	$(PRINT) "Compiling $@... "
	$(COMPILE) -o $@ -c $< -MP -MMD
	$(PRINT) "Done.\n"

### DEBUG OBJECTS ###
$(OBJECTS_D): $(BUILDPATH)/%_D.o: $(SRCPATH)/%.cpp | dirs
	$(PRINT) "Compiling $@... "
	$(DCOMPILE) -o $@ -c $< -MP -MMD
	$(PRINT) "Done.\n"

### AUTOMATIC DEPENDENCIES ###
DEPS = $(OBJECTS:%.o=%.d)
DEPS_D = $(OBJECTS_D:%.o=%.d)
-include $(DEPS)
-include $(DEPS_D)

### BINARY ###
$(BINPATH)/$(BINARY): $(OBJECTS)
	$(PRINT) "Linking $@... "
	$(COMPILE) -o $@ $(LDFLAGS) $(OBJECTS)
	$(PRINT) "Done.\n"

### DEBUG BINARY ###
$(BINPATH)/$(BINARY)_D: $(OBJECTS_D)
	$(PRINT) "Linking $@... "
	$(COMPILE) -o $@ $(LDFLAGS) $(OBJECTS_D)
	$(PRINT) "Done.\n"
