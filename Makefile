export PATH := $(shell pwd)/scripts:$(PATH)

BUILD_DIR ?= $(BUILD_DIR_BASE)_$(CONFIG)

CLANG_TIDY_ARGS = ./scripts/run-clang-tidy.py -clang-tidy-binary clang-tidy-12 -clang-apply-replacements-binary clang-apply-replacements-12 -j $(shell nproc)

.PHONY: tidy-check
tidy-check:
	${CLANG_TIDY_ARGS} -quiet -p $(BUILD_DIR) -header-filter $(shell realpath ./src) $(shell realpath ./src)

.PHONY: tidy-fix
tidy-fix:
	${CLANG_TIDY_ARGS} -fix -quiet -p $(BUILD_DIR) -header-filter $(shell realpath ./src) $(shell realpath ./src)


.PHONY: coverage
coverage:
	cmake -S . -G Ninja -B cmake-build-coverage -D CMAKE_BUILD_TYPE=Debug -D LT_BUILD_COVERAGE=TRUE
	cmake --build cmake-build-coverage
	cd cmake-build-coverage && ctest

.PHONY: coverage-html
coverage-html: coverage
	cd cmake-build-coverage && gcovr --html --html-details --exclude-unreachable-branches -o coverage.html -r ../src -j ${shell nproc} -s .

.PHONY: coverage-xml
coverage-xml: coverage
	cd cmake-build-coverage && gcovr --xml-pretty --exclude-unreachable-branches -o coverage.xml  -r ../src -j ${shell nproc} -s .

.PHONY: stats
stats:
	@cloc -by-file-by-lang --exclude-dir=3rd_party --exclude-ext=svg --vcs=git .

.PHONY: format
format:
	@find auxiliary -iname '*.hpp' -o -iname '*.h' -o -iname '*.cpp' | xargs clang-format-12 -i
	@find src -iname '*.hpp' -o -iname '*.h' -o -iname '*.cpp' | xargs clang-format-12 -i
	@find test -iname '*.hpp' -o -iname '*.h' -o -iname '*.cpp' | xargs clang-format-12 -i
	@find unitTests -iname '*.hpp' -o -iname '*.h' -o -iname '*.cpp' | xargs clang-format-12 -i

.PHONY: format-check
format-check:
	@find src -iname '*.hpp' -o -iname '*.h' -o -iname '*.cpp' | xargs -n 1 -P 1 -I{} -t sh -c 'clang-format-12 -style=file {} | diff - {}'