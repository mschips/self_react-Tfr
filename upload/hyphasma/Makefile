##############################################################################
# Hyphasma
# --------------------------------------------------------------------------- 
#
# Makefile
#
# Summary:
#   Makefile recognizes source code and generates dependencies
#   automatically:
# 
#   - each file src/*.$(SRC_EXT) is considered a source file
#   - for each source file, an build target for an object file
#     (.o) is created that depends on the source file
#   - included headers are automatically added as dependencies
#
#   Relevant folders:
#   -- ./
#     | -- ./src (.cpp and .h)
#     | -- ./build (.o)
#     | -- ./bin (binary, also linked to ./)
#
# Sebastian Binder <sebastiancbinder@gmail.com>

include defaults.mk

.PHONY: release
release: 
	$(PRINT) "Building hyphasma...\n"
	@if [ -f $(BINARY) ]; then rm $(BINARY); fi
	@$(MAKE) -C src all --no-print-directory
	@ln -sv $(BINPATH)/$(BINARY) . 
	$(PRINT) "Done.\n"

.PHONY: debug
debug:
	@$(MAKE) -C src debug --no-print-directory

.PHONY: tests
tests:
	@$(MAKE) -C tests all --no-print-directory
	@for test in $(BINPATH)/test_*; do \
		$$test; \
	done

.PHONY: clean
clean:
	@$(MAKE) -C src clean --no-print-directory
	@$(MAKE) -C tests clean --no-print-directory
	$(RM) -rfv $(BINARY)

.PHONY: update_subtrees
update_subtrees:
	git fetch googletest
	git subtree pull --prefix tests/googletest googletest master --squash

.PHONY: doc
doc:
	@doxygen
	$(MAKE) -C doc/latex
	@cp -v doc/latex/refman.pdf doc/refman.pdf

.PHONY: memcheck
memcheck:
	$(MAKE) debug
	valgrind --leak-check=yes --log-file=hyphasma.memcheck --verbose --track-origins=yes $(BINPATH)/$(BINARY)_D $(PARPATH)/leda-cellrep12

.PHONY: callgraph
callgraph:
	$(MAKE) debug
	valgrind --tool=callgrind --callgrind-out-file=hyphasma.callgrind $(BINPATH)/$(BINARY)_D $(PARPATH)/leda-cellrep12
	kcachegrind hyphasma.callgrind

.PHONY: analyze
analyze: _clang-analyze
	@mv -v $(shell ls -d analyze/*/)* analyze
	@scan-view analyze

.PHONY: _clang-analyze
_clang-analyze:
	@$(MAKE) clean
	@rm -rf analyze
	@mkdir -pv analyze
	@scan-build -o analyze $(MAKE) debug -j4
