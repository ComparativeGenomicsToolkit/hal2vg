# Use together with the existing build configuration; does not rebuild hal2vg.
# make -f Makefile -f verification.mk hal2vg-fingerprint

# OpenSSL, for the digests.  Homebrew's is keg-only, so on Mac it is not on the default paths.
cryptoFlags =
ifeq (${SYS},Darwin)
	opensslPrefix := $(shell brew --prefix openssl 2>/dev/null)
	ifneq (${opensslPrefix},)
		cryptoFlags = -I${opensslPrefix}/include -L${opensslPrefix}/lib
	endif
endif

hal2vg-fingerprint: tests/fingerprint_graph.cpp ${basicLibsDependencies}
	${cpp} ${CXXFLAGS} ${cryptoFlags} ${ompLibs} -pthread -I . tests/fingerprint_graph.cpp ${basicLibs} -lcrypto -o hal2vg-fingerprint

fingerprint-fixture: tests/fingerprint_fixture.cpp ${basicLibsDependencies}
	${cpp} ${CXXFLAGS} ${ompLibs} -pthread -I . tests/fingerprint_fixture.cpp ${basicLibs} -o fingerprint-fixture

.PHONY: test-graph-fingerprints
test-graph-fingerprints: hal2vg-fingerprint fingerprint-fixture
	python3 tests/test_graph_fingerprints.py
